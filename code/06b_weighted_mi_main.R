# 06b_weighted_mi_main.R -- final numbers for Table 2 / Table 3 under the weighted-MI estimator
# Produces, for both cut-off definitions:
#   (a) age-sex adjusted RRRs   (b) fully adjusted RRRs   (c) ordinal trend per low reading
#   (d) standardized absolute risks at 60-80 y
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")

mn_vcov <- function(beta, X, y, w, cluster){
  K <- 2; p <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P/rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y==k))
  A <- matrix(0, K*p, K*p)
  Wm <- list(P[,2]*(1-P[,2]), P[,3]*(1-P[,3]), -P[,2]*P[,3])
  bl <- matrix(c(1,3,3,2), 2, 2, byrow=TRUE)
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[bl[k,l]]]
    A[((k-1)*p+1):(k*p), ((l-1)*p+1):(l*p)] <- crossprod(X, X*ww)
  }
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[, ((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1])) * X
  Sc <- rowsum(S, group=cluster, reorder=FALSE); G <- nrow(Sc)
  Ai <- tryCatch(solve(A), error=function(e){ cat("  [warn] singular A, pseudo-inverse used\n"); MASS::ginv(A) })
  list(V = (G/(G-1))*((N-1)/(N-K*p)) * Ai %*% crossprod(Sc) %*% Ai)
}

mkX <- function(d, covars, expo){
  X <- cbind(1, as.matrix(d[, covars]))
  if (expo == "grp") for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  if (expo == "nr")  X <- cbind(X, d$nr_low)
  colnames(X) <- c("const", covars,
    if (expo=="grp") paste0("grp", c("remitted","incident","persistent")) else "nr_low")
  X
}

rubin <- function(bs, vs){
  mm <- length(bs); b <- Reduce(`+`, bs)/mm
  W <- Reduce(`+`, vs)/mm
  B <- Reduce(`+`, lapply(bs, function(z) (z-b)^2))/(mm-1)
  T <- W + (1+1/mm)*B
  lam <- (1+1/mm)*B/pmax(T,1e-12)
  df <- (mm-1)/pmax(lam,1e-12)^2
  list(b=b, se=sqrt(T), df=ifelse(is.finite(df) & df>0, df, 1e5))
}

pool_model <- function(d0, imp, covars, expo){
  bs <- list(); vs <- list()
  for (i in 1:20) {
    di <- complete(imp, i)
    di$y <- factor(di$y, levels=c("alive_free","disabled","dead"))
    di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
    di$nr_low <- di$nr
    X <- mkX(di, covars, expo)
    ex <- if (expo=="grp") "+grp" else "+nr_low"
    f <- as.formula(paste("y ~", paste(c(covars, if(expo=="grp") "grp" else "nr_low"), collapse="+")))
    fit <- multinom(f, data=di, weights=di$wsc, trace=FALSE, maxit=500)
    bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- colnames(X)
    vc <- mn_vcov(bb, X, as.integer(di$y)-1, di$wsc, di$communityID)
    bs[[i]] <- as.numeric(t(bb)); vs[[i]] <- as.numeric(diag(vc$V))
  }
  pr <- rubin(bs, vs); p <- ncol(X)
  B <- matrix(pr$b, 2, byrow=TRUE); S <- matrix(pr$se, 2, byrow=TRUE); D <- matrix(pr$df, 2, byrow=TRUE)
  colnames(B) <- colnames(S) <- colnames(X)
  list(B=B, S=S, D=D, cols=colnames(X), p=p)
}

report <- function(res, label){
  B <- res$B; S <- res$S; D <- res$D
  out <- list()
  for (k in 1:2) for (j in (res$p-2):res$p) {
    nm <- res$cols[j]
    g <- sub("^grp","",nm)
    b <- B[k,j]; s <- S[k,j]; df <- D[k,j]
    out[[length(out)+1]] <- data.frame(model=label, outcome=c("disabled","dead")[k], group=g,
      RRR=exp(b), lo=exp(b-qt(.975,df)*s), hi=exp(b+qt(.975,df)*s), p=2*pt(-abs(b/s),df),
      stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}

all_rrr <- list(); all_abs <- list()
for (cut in c("19","25")) {
  b0 <- get_base(m, paste0("grp",cut))
  d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13>0) %>%
        mutate(nr = low11_19 + low13_19,
               nr_low = ifelse(cut=="19", low11_19 + low13_19, low11_25 + low13_25),
               edu_f = factor(edu), wsc = w13/mean(w13, na.rm=TRUE))
  d0$nr_low <- d0$nr_low - 0
  cat(sprintf("\n### cut-offs %s | n=%d | clusters=%d\n", cut, nrow(d0), length(unique(d0$communityID))))
  imp <- mice(d0[, c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y","nr_low")],
              m=20, seed=20261023, printFlag=FALSE)
  pool1 <- function(covars, expo){
    bs <- list(); vs <- list()
    for (i in 1:20) {
      di <- complete(imp, i)
      di$y <- factor(di$y, levels=c("alive_free","disabled","dead"))
      di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
      di$edu <- as.numeric(as.character(di$edu_f))
      for (kk in 1:4) di[[paste0("edu",kk)]] <- as.numeric(di$edu==kk)
      di$communityID <- d0$communityID; di$wsc <- d0$wsc
      X <- mkX(di, covars, expo)
      f <- as.formula(paste("y ~", paste(c(covars, if (expo=="grp") "grp" else "nr_low"), collapse="+")))
      fit <- multinom(f, data=di, weights=di$wsc, trace=FALSE, maxit=500)
      bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- colnames(X)
      vc <- mn_vcov(bb, X, as.integer(di$y)-1, di$wsc, di$communityID)
      bs[[i]] <- as.numeric(t(bb)); vs[[i]] <- as.numeric(diag(vc$V))
    }
    pr <- rubin(bs, vs); p <- ncol(X)
    list(B=matrix(pr$b,2,byrow=TRUE), S=matrix(pr$se,2,byrow=TRUE),
         D=matrix(pr$df,2,byrow=TRUE), cols=colnames(X), p=p)
  }
  rA <- pool1(c("age10","male"), "grp");       a <- report(rA, paste0("age+sex,", cut))
  rF <- pool1(COV, "grp");                     b <- report(rF, paste0("fully adj,", cut))
  rT <- pool1(COV, "nr");                      tr <- c(RRR=exp(rT$B[1,rT$p]), lo=exp(rT$B[1,rT$p]-qt(.975,rT$D[1,rT$p])*rT$S[1,rT$p]),
                                                     hi=exp(rT$B[1,rT$p]+qt(.975,rT$D[1,rT$p])*rT$S[1,rT$p]),
                                                     RRR_x=exp(rT$B[2,rT$p]),
                                                     lo_x=exp(rT$B[2,rT$p]-qt(.975,rT$D[2,rT$p])*rT$S[2,rT$p]),
                                                     hi_x=exp(rT$B[2,rT$p]+qt(.975,rT$D[2,rT$p])*rT$S[2,rT$p]))
  cat(sprintf("ordinal per additional low reading (%s): disability %.2f (%.2f-%.2f); death %.2f (%.2f-%.2f)\n",
      cut, tr[1], tr[2], tr[3], tr[4], tr[5], tr[6]))
  print(rbind(a, b), row.names=FALSE, digits=2)
  all_rrr[[cut]] <- rbind(a, b)
  # absolute risks
  prof <- c(1, (mean(d0$age13)-65)/10, 0.5, mean(d0$bmi,na.rm=TRUE), mean(d0$nchronic,na.rm=TRUE),
            mean(d0$cesd,na.rm=TRUE), mean(d0$smoke,na.rm=TRUE), 0,0,0,0)
  names(prof) <- rF$cols[1:11]
  ab <- list()
  for (age in c(60,65,70,75,80)) for (g in GRP) {
    x <- prof; x["age10"] <- (age-65)/10
    x[paste0("grp", c("remitted","incident","persistent"))] <- 0
    if (g!="never") x[paste0("grp",g)] <- 1
    e <- as.numeric(rF$B %*% x); pp <- exp(c(0,e)); pp <- pp/sum(pp)
    ab[[length(ab)+1]] <- data.frame(cutoff=cut, age=age, group=g, disabled=100*pp[2], dead=100*pp[3])
  }
  all_abs[[cut]] <- do.call(rbind, ab)
  cat("absolute risk:\n"); print(reshape(all_abs[[cut]][,c("age","group","disabled","dead")],
      idvar="age", timevar="group", direction="wide"), row.names=FALSE, digits=1)
}
rrr <- do.call(rbind, all_rrr); abs_r <- do.call(rbind, all_abs)
write.csv(rrr, "out/S2_weighted_MI_rrr.csv", row.names=FALSE)
write.csv(abs_r, "out/S2_weighted_MI_absolute_risk.csv", row.names=FALSE)
cat("\nwritten out/S2_weighted_MI_rrr.csv, S2_weighted_MI_absolute_risk.csv\n")
