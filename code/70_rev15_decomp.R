# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 70_rev15_decomp.R  (v2: with and without adjustment for the 2013 grip value)
# P1-3: decomposition of the persistent-vs-never-low contrast
#          estimand (multinomial RRR vs binary OR)
#        x sample   (multiply imputed vs complete case)
#        x weighting(design weight vs unweighted)
#        x adjustment (with / without 2013 grip per 5 kg)
# Fidelity checks: primary must be 1.83 / 2.89; adjusted binary CC unweighted must be ~1.18 / 1.28 (S3 B2).
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr); library(sandwich)})
set.seed(20261023)
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
    A[((k-1)*p+1):(k*p), ((l-1)*p+1):((l-1)*p+p)] <- crossprod(X, X*ww)
  }
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[, ((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1])) * X
  Sc <- rowsum(S, group=cluster, reorder=FALSE); G <- nrow(Sc)
  Ai <- tryCatch(solve(A), error=function(e) MASS::ginv(A))
  list(V = (G/(G-1))*((N-1)/(N-K*p)) * Ai %*% crossprod(Sc) %*% Ai)
}

mkX <- function(d, covars){
  X <- cbind(1, as.matrix(d[, covars]))
  for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  colnames(X) <- c("const", covars, paste0("grp", c("remitted","incident","persistent")))
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

addEdu <- function(d){
  d$edu <- as.numeric(as.character(d$edu_f))
  for (kk in 1:4) d[[paste0("edu",kk)]] <- as.numeric(d$edu==kk)
  d
}

b0 <- get_base(m, "grp19")
d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13>0) %>%
      mutate(wsc = w13/mean(w13, na.rm=TRUE), edu_f = factor(edu), g13_5 = g13/5)
d0$grp <- relevel(factor(d0$grp, levels=GRP), ref="never")
d0$y   <- factor(d0$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
cat(sprintf("analysis set n = %d | clusters = %d | deaths = %d | disabled = %d\n",
            nrow(d0), length(unique(d0$communityID)), sum(d0$y=="dead"), sum(d0$y=="disabled")))

IMP <- mice(d0[, c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y","g13_5")],
            m=20, seed=20261023, printFlag=FALSE)
cat("imputations done\n")

res <- list()
push <- function(est, lo, hi, estlab, sample, weight, adj, out){
  res[[length(res)+1]] <<- data.frame(estimand=estlab, sample=sample, weighting=weight,
    adjustment=ifelse(adj,"+ 2013 grip","unadjusted"), outcome=out,
    est=est, lo=lo, hi=hi, stringsAsFactors=FALSE)
}
FO <- function(adj) paste("y ~ age10+male+bmi+nchronic+cesd+smoke+edu1+edu2+edu3+edu4",
                          if (adj) "+g13_5" else "", "+grp")
FB <- function(adj) paste("yb ~ age10+male+bmi+nchronic+cesd+smoke+edu1+edu2+edu3+edu4",
                          if (adj) "+g13_5" else "", "+grp")
CV <- function(adj) c(COV, if (adj) "g13_5")

mults_cc <- function(weighted, adj){
  dat <- addEdu(d0 %>% filter(!is.na(age10),!is.na(male),!is.na(bmi),!is.na(nchronic),
                              !is.na(cesd),!is.na(smoke),!is.na(edu_f),!is.na(g13)))
  dat$w <- if (weighted) dat$wsc else rep(1, nrow(dat))
  X <- mkX(dat, CV(adj))
  fit <- multinom(as.formula(FO(adj)), data=dat, weights=dat$w, trace=FALSE, maxit=500)
  bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- colnames(X)
  vc <- mn_vcov(bb, X, as.integer(dat$y)-1, dat$w, dat$communityID)
  j <- which(colnames(X)=="grppersistent")
  for (k in 1:2) {
    est <- exp(bb[k,j]); se <- sqrt(vc$V[(k-1)*ncol(X)+j, (k-1)*ncol(X)+j])
    push(est, exp(log(est)-1.96*se), exp(log(est)+1.96*se), "multinomial RRR",
         "complete case", ifelse(weighted,"weighted","unweighted"), adj, c("disability","death")[k])
  }
}
bin_cc <- function(weighted, adj){
  dat <- addEdu(d0 %>% filter(!is.na(age10),!is.na(male),!is.na(bmi),!is.na(nchronic),
                              !is.na(cesd),!is.na(smoke),!is.na(edu_f),!is.na(g13)))
  dat$w <- if (weighted) dat$wsc else rep(1, nrow(dat))
  for (k in 1:2) {
    keep <- as.integer(dat$y) %in% c(1, 1+k)
    dk <- dat[keep, ]; dk$yb <- as.numeric(as.integer(dk$y)==k+1)
    g <- suppressWarnings(glm(as.formula(FB(adj)), data=dk, weights=dk$w, family=binomial()))
    vv <- tryCatch(sandwich::vcovCL(g, cluster=dk$communityID, type="HC0"), error=function(e) vcov(g))
    est <- coef(g)["grppersistent"]; se <- sqrt(vv["grppersistent","grppersistent"])
    push(exp(est), exp(est-1.96*se), exp(est+1.96*se), "binary OR", "complete case",
         ifelse(weighted,"weighted","unweighted"), adj, c("disability","death")[k])
  }
}
for (adj in c(FALSE, TRUE)) { mults_cc(FALSE, adj); mults_cc(TRUE, adj); bin_cc(FALSE, adj); bin_cc(TRUE, adj) }

pool_imp <- function(weighted, adj){
  bs <- list(); vs <- list(); bn <- list(); vn <- list()
  for (i in 1:20) {
    di <- addEdu(complete(IMP, i))
    di$communityID <- d0$communityID; di$wsc <- d0$wsc
    di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
    di$w <- if (weighted) d0$wsc else rep(1, nrow(di))
    X <- mkX(di, CV(adj))
    fit <- multinom(as.formula(FO(adj)), data=di, weights=di$w, trace=FALSE, maxit=500)
    bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- colnames(X)
    vc <- mn_vcov(bb, X, as.integer(di$y)-1, di$w, di$communityID)
    j <- which(colnames(X)=="grppersistent")
    bs[[i]] <- bb[,j]
    vs[[i]] <- c(vc$V[j,j], vc$V[ncol(X)+j, ncol(X)+j])
    for (k in 1:2) {
      keep <- as.integer(di$y) %in% c(1, 1+k)
      dk <- di[keep, ]; dk$yb <- as.numeric(as.integer(dk$y)==k+1)
      g <- suppressWarnings(glm(as.formula(FB(adj)), data=dk, weights=dk$w, family=binomial()))
      vv <- tryCatch(sandwich::vcovCL(g, cluster=dk$communityID, type="HC0"), error=function(e) vcov(g))
      bn[[length(bn)+1]] <- coef(g)["grppersistent"]
      vn[[length(vn)+1]] <- vv["grppersistent","grppersistent"]
    }
  }
  pm <- rubin(bs, vs)
  for (k in 1:2) push(exp(pm$b[k]), exp(pm$b[k]-qt(.975,pm$df[k])*pm$se[k]),
                      exp(pm$b[k]+qt(.975,pm$df[k])*pm$se[k]), "multinomial RRR", "imputed",
                      ifelse(weighted,"weighted","unweighted"), adj, c("disability","death")[k])
  for (k in 1:2) {
    pb <- rubin(lapply(seq(k, length(bn), by=2), function(z) bn[[z]]),
                lapply(seq(k, length(vn), by=2), function(z) vn[[z]]))
    push(exp(pb$b), exp(pb$b-qt(.975,pb$df)*pb$se), exp(pb$b+qt(.975,pb$df)*pb$se), "binary OR",
         "imputed", ifelse(weighted,"weighted","unweighted"), adj, c("disability","death")[k])
  }
}
for (adj in c(FALSE, TRUE)) { pool_imp(TRUE, adj); pool_imp(FALSE, adj) }

out <- do.call(rbind, res)
out$cell <- sprintf("%.2f (%.2f-%.2f)", out$est, out$lo, out$hi)
out <- out[order(out$outcome, out$adjustment, out$estimand, out$sample, out$weighting), ]
print(out[, c("outcome","adjustment","estimand","sample","weighting","cell")], row.names=FALSE)
write.csv(out, "out/70_decomp.csv", row.names=FALSE)

sink("out/70_decomp.txt")
cat("P1-3 decomposition | persistent vs never-low\n")
cat(sprintf("n = %d; clusters = %d; deaths = %d; ADL disability = %d\n\n",
            nrow(d0), length(unique(d0$communityID)), sum(d0$y=="dead"), sum(d0$y=="disabled")))
print(out[, c("outcome","adjustment","estimand","sample","weighting","cell")], row.names=FALSE)
cat("\n--- marginal effect of each feature, within outcome x adjustment ---\n")
for (oc in c("disability","death")) for (ad in unique(out$adjustment)) {
  dd <- out[out$outcome==oc & out$adjustment==ad, ]
  if (!nrow(dd)) next
  em <- tapply(dd$est, dd$estimand, mean); es <- tapply(dd$est, dd$sample, mean); ew <- tapply(dd$est, dd$weighting, mean)
  cat(sprintf("%-10s %-11s | estimand %s (ratio %.2f) | sample %s (ratio %.2f) | weighting %s (ratio %.2f)\n",
      oc, ad,
      paste(sprintf("%.2f", em), collapse=" vs "), max(em)/min(em),
      paste(sprintf("%.2f", es), collapse=" vs "), max(es)/min(es),
      paste(sprintf("%.2f", ew), collapse=" vs "), max(ew)/min(ew)))
}
sink()
cat("\nwritten out/70_decomp.csv / .txt\n")
