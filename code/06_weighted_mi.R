# 06_weighted_mi.R -- weighted multiple-imputation estimates for the main table
# Keeps the manuscript's MI approach (m = 20) but adds CHARLS 2013 individual weights
# and community-clustered variance, so Table 2 and Table 3 share one estimator.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")

CU <- function(cut) if (cut=="19") c(28,18) else NULL

mn_cluster_vcov <- function(beta, X, y, w, cluster){
  K <- 2; p <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P/rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y==k))
  A <- matrix(0, K*p, K*p)
  Wm <- list(P[,2]*(1-P[,2]), P[,3]*(1-P[,3]), -P[,2]*P[,3])
  blocks <- matrix(c(1,3,3,2), 2, 2, byrow=TRUE)
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[blocks[k,l]]]
    A[((k-1)*p+1):(k*p), ((l-1)*p+1):(l*p)] <- crossprod(X, X*ww)
  }
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[, ((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1])) * X
  Sc <- rowsum(S, group=cluster, reorder=FALSE)
  G <- nrow(Sc); Ainv <- solve(A)
  cfac <- G/(G-1) * (N-1)/(N-K*p)
  list(V = cfac * Ainv %*% crossprod(Sc) %*% Ainv, G=G, p=K*p)
}

build_X <- function(d, covars){
  X <- cbind(1, as.matrix(d[, covars]))
  for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  colnames(X) <- c("const", covars, paste0("grp", c("remitted","incident","persistent")))
  X
}

rubin <- function(betas, vars, dfcom){
  m <- length(betas); b <- Reduce(`+`, betas)/m
  W <- Reduce(`+`, vars)/m
  B <- Reduce(`+`, lapply(betas, function(z) (z-b)^2))/(m-1)
  T <- W + (1+1/m)*B
  lam <- (1+1/m)*B/T
  df <- (m-1)/pmax(lam, 1e-10)^2
  df <- ifelse(is.finite(df) & df>0, df, dfcom)
  list(b=b, se=sqrt(T), df=df, lam=lam)
}

pool_cut <- function(cut, label){
  b0 <- get_base(m, paste0("grp",cut))
  d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13>0) %>%
        mutate(g13_5 = g13/5, edu_f = factor(edu), wsc = w13/mean(w13, na.rm=TRUE))
  cat(sprintf("\n=== cut-offs %s (%s): n = %d, deaths = %d, clusters = %d ===\n",
      cut, label, nrow(d0), sum(d0$oc==2), length(unique(d0$communityID))))
  vars <- c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y")
  imp_dat <- d0[, c(vars, "communityID","wsc")]
  imp <- mice(imp_dat, m=20, seed=20261023, printFlag=FALSE, remove.data=TRUE)
  betas <- list(); vars_list <- list()
  for (i in 1:20) {
    di <- complete(imp, i)
    di$y <- factor(di$y, levels=c("alive_free","disabled","dead"))
    di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
    di$edu <- as.numeric(as.character(di$edu_f))
    for (k in 1:4) di[[paste0("edu",k)]] <- as.numeric(di$edu==k)
    X <- build_X(di, COV)
    f <- as.formula(paste("y ~", paste(c(COV,"grp"), collapse="+")))
    fit <- multinom(f, data=di, weights=di$wsc, trace=FALSE, maxit=500)
    bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- colnames(X)
    vc <- mn_cluster_vcov(bb, X, as.integer(di$y)-1, di$wsc, di$communityID)
    betas[[i]] <- as.numeric(t(bb))          # row-major: eq1 all coefs, then eq2
    vars_list[[i]] <- as.numeric(diag(vc$V))
  }
  p <- ncol(X)
  pr <- rubin(betas, vars_list, dfcom=nrow(d0)-2*p)
  B <- matrix(pr$b, nrow=2, byrow=TRUE); SE <- matrix(pr$se, nrow=2, byrow=TRUE)
  DF <- matrix(pr$df, nrow=2, byrow=TRUE)
  colnames(B) <- colnames(SE) <- colnames(X)
  rows <- list()
  for (k in 1:2) for (g in c("remitted","incident","persistent")) {
    j <- which(colnames(X)==paste0("grp",g)); b <- B[k,j]; s <- SE[k,j]; df <- DF[k,j]
    rows[[length(rows)+1]] <- data.frame(cutoff=cut, outcome=c("disabled","dead")[k], group=g,
      RRR=exp(b), lo=exp(b-qt(0.975, df)*s), hi=exp(b+qt(0.975, df)*s),
      p=2*pt(-abs(b/s), df), n=nrow(d0), clusters=length(unique(d0$communityID)),
      stringsAsFactors=FALSE)
  }
  rows <- do.call(rbind, rows)
  cat(sprintf("cluster-robust G = %d | %s\n", length(unique(d0$communityID)),
      paste(sprintf("%s %s %.2f (%.2f-%.2f)", rows$outcome, rows$group, rows$RRR, rows$lo, rows$hi), collapse=" | ")))
  # absolute risks from pooled betas, weighted mean covariate profile
  prof <- c(1, (mean(d0$age13)-65)/10, mean(d0$male, na.rm=TRUE), mean(d0$bmi, na.rm=TRUE),
            mean(d0$nchronic, na.rm=TRUE), mean(d0$cesd, na.rm=TRUE),
            mean(d0$smoke, na.rm=TRUE), 0, 0, 0, 0)
  names(prof) <- colnames(X)[1:11]
  abs_r <- list()
  for (a in c(60,65,70,75,80)) for (g in GRP) {
    x <- prof; x["age10"] <- (a-65)/10; x["male"] <- 0.5
    x[paste0("grp", c("remitted","incident","persistent"))] <- 0
    if (g!="never") x[paste0("grp",g)] <- 1
    eta <- as.numeric(B %*% x); pp <- exp(c(0, eta)); pp <- pp/sum(pp)
    abs_r[[length(abs_r)+1]] <- data.frame(cutoff=cut, age=a, group=g,
      disabled=100*pp[2], dead=100*pp[3], stringsAsFactors=FALSE)
  }
  list(rrr=rows, abs=do.call(rbind, abs_r), imp=imp)
}

r19 <- pool_cut("19", "AWGS 2019, primary")
r25 <- pool_cut("25", "AWGS 2025, sensitivity")
rrr <- rbind(r19$rrr, r25$rrr); abs_r <- rbind(r19$abs, r25$abs)
write.csv(rrr, "out/S2_weighted_MI_rrr.csv", row.names=FALSE)
write.csv(abs_r, "out/S2_weighted_MI_absolute_risk.csv", row.names=FALSE)

cat("\n--- absolute risk, weighted MI model (%) ---\n")
for (k in c("19","25")) {
  cat("cut-offs", k, "\n")
  w <- reshape(abs_r[abs_r$cutoff==k, c("age","group","disabled","dead")],
               idvar="age", timevar="group", direction="wide")
  print(w, row.names=FALSE, digits=1)
}
saveRDS(list(rrr=rrr, abs=abs_r), "out/S2_weighted_MI.rds")
cat("\nDone -> out/S2_weighted_MI_rrr.csv, out/S2_weighted_MI_absolute_risk.csv\n")
