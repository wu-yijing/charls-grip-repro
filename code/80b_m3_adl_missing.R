# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# M3: the participants whose 2013 ADL block was entirely missing. Variant A excludes them
# (the primary analysis); variant B counts them as ADL-free. Weighted, multiply imputed,
# cluster-robust multinomial models. Produces Supplementary Table S23.
# 80b_m3_adl_missing.R -- M3: the 2,617 with an entirely missing 2013 ADL block
# Variant A: primary (missing ADL -> excluded)
# Variant B: missing ADL -> ADL-free (n should rise to ~8,669)
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
  Ai <- tryCatch(solve(A), error=function(e) MASS::ginv(A))
  list(V = (G/(G-1))*((N-1)/(N-K*p)) * Ai %*% crossprod(Sc) %*% Ai)
}
rubin <- function(bs, vs){
  mm <- length(bs); b <- Reduce(`+`, bs)/mm; W <- Reduce(`+`, vs)/mm
  B <- Reduce(`+`, lapply(bs, function(z) (z-b)^2))/(mm-1)
  T <- W + (1+1/mm)*B; lam <- (1+1/mm)*B/pmax(T,1e-12); df <- (mm-1)/pmax(lam,1e-12)^2
  list(b=b, se=sqrt(T), df=ifelse(is.finite(df) & df>0, df, 1e5))
}
mkX <- function(d, covars){
  X <- cbind(1, as.matrix(d[, covars]))
  for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  colnames(X) <- c("const", covars, paste0("grp", c("remitted","incident","persistent")))
  X
}

run_cohort <- function(include_missing, label){
  b <- m %>% filter(if (include_missing) (a13 %in% 0 | is.na(a13)) else (a13 %in% 0),
                    !is.na(grp19), !is.na(age13), !is.na(male))
  b$grp <- relevel(factor(b$grp19, levels=GRP), ref="never")
  b$y   <- factor(b$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
  d0 <- b %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
        mutate(nr_low = low11_19 + low13_19, edu_f = factor(edu))
  d0$wsc <- d0$w13/mean(d0$w13, na.rm=TRUE)
  cat(sprintf("\n### %s | eligible=%d | modelled=%d | deaths=%d | disabled=%d | clusters=%d\n",
      label, nrow(b), nrow(d0), sum(d0$y=="dead"), sum(d0$y=="disabled"),
      length(unique(d0$communityID))))
  cat("    group sizes (modelled): "); print(table(d0$grp))
  imp <- mice(d0[, c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y","nr_low")],
              m=20, seed=20261023, printFlag=FALSE)
  for (spec in list(list(cv=c("age10","male"), nm="age+sex"), list(cv=COV, nm="fully adjusted"))) {
    covars <- spec$cv; bs <- list(); vs <- list(); Xcols <- NULL
    for (i in 1:20) {
      di <- complete(imp, i)
      di$y <- factor(di$y, levels=c("alive_free","disabled","dead"))
      di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
      di$edu <- as.numeric(as.character(di$edu_f))
      for (kk in 1:4) di[[paste0("edu",kk)]] <- as.numeric(di$edu==kk)
      di$communityID <- d0$communityID; di$wsc <- d0$wsc
      X <- mkX(di, covars); Xcols <- colnames(X)
      f <- as.formula(paste("y ~", paste(c(covars,"grp"), collapse="+")))
      fit <- multinom(f, data=di, weights=di$wsc, trace=FALSE, maxit=500)
      bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- Xcols
      vc <- mn_vcov(bb, X, as.integer(di$y)-1, di$wsc, di$communityID)
      bs[[i]] <- as.numeric(t(bb)); vs[[i]] <- as.numeric(diag(vc$V))
    }
    pr <- rubin(bs, vs); p <- length(Xcols)
    B <- matrix(pr$b,2,byrow=TRUE); S <- matrix(pr$se,2,byrow=TRUE); D <- matrix(pr$df,2,byrow=TRUE)
    for (k in 1:2) for (j in (p-2):p) {
      nm <- Xcols[j]; dfk <- D[k,j]
      cat(sprintf("  %-15s %-8s %-10s %.2f (%.2f-%.2f)  p=%.3f\n",
        spec$nm, c("disability","death")[k], sub("^grp","",nm),
        exp(B[k,j]), exp(B[k,j]-qt(.975,dfk)*S[k,j]), exp(B[k,j]+qt(.975,dfk)*S[k,j]),
        2*pt(-abs(B[k,j]/S[k,j]), dfk)))
    }
  }
}

run_cohort(FALSE, "A: primary (missing ADL excluded)")
run_cohort(TRUE,  "B: missing ADL counted as ADL-free")
