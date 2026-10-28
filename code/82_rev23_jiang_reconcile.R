# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
#
# Axis-by-axis reconciliation of the persistent-versus-never mortality contrast with
# Jiang 2026 (reference 17): cut-off set, age restriction, estimand (Section 4.2).
# 82_rev23_jiang_reconcile.R -- align our contrast with Jiang 2026 (ref 17) axis by axis
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages({library(dplyr); library(nnet)})
m <- readRDS("out/jamda_full.rds")

mn_vcov <- function(beta, X, y, w, cluster){
  K <- 2; p <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P/rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y==k))
  A <- matrix(0, K*p, K*p); Wm <- list(P[,2]*(1-P[,2]), P[,3]*(1-P[,3]), -P[,2]*P[,3])
  bl <- matrix(c(1,3,3,2),2,2,byrow=TRUE)
  for (k in 1:K) for (l in 1:K) { ww <- w*Wm[[bl[k,l]]]; A[((k-1)*p+1):(k*p),((l-1)*p+1):(l*p)] <- crossprod(X, X*ww) }
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[,((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1]))*X
  Sc <- rowsum(S, group=cluster, reorder=FALSE); G <- nrow(Sc)
  Ai <- tryCatch(solve(A), error=function(e) MASS::ginv(A))
  list(V=(G/(G-1))*((N-1)/(N-K*p))*Ai%*%crossprod(Sc)%*%Ai)
}
prep <- function(col, minage=0){
  b <- m %>% filter(a13 %in% 0, !is.na(.data[[col]]), !is.na(age13), !is.na(male), age13>=minage)
  b$grp <- relevel(factor(b[[col]], levels=GRP), ref="never")
  b$y <- factor(b$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
  b$wsc <- b$w13/mean(b$w13, na.rm=TRUE); b
}
fit1 <- function(b, covars=COV){
  d <- b[complete.cases(b[, c("y",covars,"grp","wsc","communityID")]), ]
  f <- as.formula(paste("y ~", paste(c(covars,"grp"),collapse="+")))
  fit <- multinom(f, data=d, weights=d$wsc, trace=FALSE, maxit=500)
  X <- cbind(1, as.matrix(d[,covars]))
  for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  colnames(X) <- c("const",covars,paste0("grp",c("remitted","incident","persistent")))
  bb <- matrix(coef(fit),2,byrow=FALSE); V <- mn_vcov(bb,X,as.integer(d$y)-1,d$wsc,d$communityID)$V
  j <- which(colnames(X)=="grppersistent"); p <- ncol(X)
  est <- bb[2,j]; se <- sqrt(V[(p)+j,(p)+j])
  data.frame(n=nrow(d), deaths=sum(d$y=="dead"), RRR=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se))
}
rows <- list(
  c("AWGS 2019 cut-offs, all ages (primary)", "19", "0"),
  c("AWGS 2025 cut-offs, all ages", "25", "0"),
  c("AWGS 2019 cut-offs, age >= 60", "19", "60"),
  c("AWGS 2025 cut-offs, age >= 60 (aligned with ref 17)", "25", "60"),
  c("AWGS 2019 cut-offs, age >= 65", "19", "65"))
out <- do.call(rbind, lapply(rows, function(r){
  x <- fit1(prep(paste0("grp", r[2]), as.numeric(r[3]))); cbind(spec=r[1], x)
}))
print(out, row.names=FALSE, digits=3)
write.csv(out, "out/rev23_jiang_reconcile.csv", row.names=FALSE)
