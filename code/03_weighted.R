# 03_weighted.R -- Analysis (2): CHARLS sampling weights + stratified cluster design
# Reviewer M3 / M5. Re-estimates the SAME multinomial logit with:
#   (a) official individual weights (2013 cross-sectional, non-response adjusted)
#   (b) community (=PSU) cluster-robust variance
#   (c) full design-based (strata x PSU x weights) binomial contrasts for comparison
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(survey)})
options(survey.lonely.psu = "adjust")   # some strata contain a single community
m <- readRDS("out/jamda_full.rds")
sink("out/S2_weighted_report.txt", split=TRUE)
cat("==================================================================\n")
cat("SUPPLEMENTARY TABLE S2 (generated", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
cat("Analysis 2 (M3/M5): sampling weights and design-based variance\n")
cat("==================================================================\n\n")

## ---- cluster-robust variance for a multinomial logit (3 categories, ref = alive-free)
mn_cluster_vcov <- function(fit, X, y, w, cluster){
  K <- 2; p <- ncol(X); N <- nrow(X)
  B <- matrix(coef(fit), nrow=K, ncol=p, byrow=FALSE)   # K x p
  eta <- X %*% t(B)
  P <- exp(cbind(0, eta)); P <- P/rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y==k+1))     # disabled, dead
  # bread: A = sum_i w_i * W_i (x) x x'
  A <- matrix(0, K*p, K*p)
  Wm <- list()
  Wm[[1]] <- P[,2]*(1-P[,2]); Wm[[2]] <- P[,3]*(1-P[,3]); Wm[[3]] <- -P[,2]*P[,3]
  blocks <- matrix(0, K, K); blocks[1,1] <- 1; blocks[2,2] <- 2; blocks[1,2] <- 3; blocks[2,1] <- 3
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[blocks[k,l]]]
    A[((k-1)*p+1):(k*p), ((l-1)*p+1):(l*p)] <- crossprod(X, X*ww)
  }
  # scores
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[, ((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1])) * X
  Sc <- rowsum(S, group=cluster, reorder=FALSE)
  G <- nrow(Sc)
  Bm <- crossprod(Sc)
  Ainv <- solve(A)
  cfac <- G/(G-1) * (N-1)/(N-K*p)
  V <- cfac * Ainv %*% Bm %*% Ainv
  list(V=V, G=G, N=N, beta=B)
}

build_X <- function(d, covars){
  X <- cbind(1, as.matrix(d[, covars]))
  for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(d$grp==g))
  colnames(X) <- c("const", covars, paste0("grp", c("remitted","incident","persistent")))
  X
}

mk_tab <- function(beta, se, label, n, extra=""){
  rows <- list()
  for (k in 1:2) {
    lab <- c("disabled","dead")[k]
    for (g in c("remitted","incident","persistent")) {
      j <- which(colnames(beta)==paste0("grp",g))
      b <- beta[k,j]; s <- se[k,j]
      rows[[length(rows)+1]] <- data.frame(model=label, outcome=lab, group=g,
        RRR=exp(b), lo=exp(b-1.96*s), hi=exp(b+1.96*s), p=2*pnorm(-abs(b/s)),
        n=n, note=extra, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

res_compare <- list()
for (cut in c("19","25")) {
  cat("###### cut-offs:", ifelse(cut=="19","AWGS 2019 (submitted)","AWGS 2025"), "######\n\n")
  b <- get_base(m, paste0("grp",cut))
  # weight availability
  cat("cohort n =", nrow(b), "| w13 available:", sum(!is.na(b$w13) & b$w13>0),
      "| w11:", sum(!is.na(b$w11) & b$w11>0),
      "| w13_long:", sum(!is.na(b$w13_long) & b$w13_long>0),
      "| communities:", length(unique(na.omit(b$communityID))),
      "| strata:", length(unique(na.omit(b$stratum))), "\n")
  d <- b[complete.cases(b[, c("y", COV, "grp")]) & !is.na(b$w13) & b$w13>0, ]
  cat("complete-case with weight: n =", nrow(d), "| clusters:",
      length(unique(d$communityID)), "\n")
  X <- build_X(d, COV)
  w <- d$w13; w <- w/mean(w)   # rescale for numerical stability; ratios unchanged
  f <- as.formula(paste("y ~", paste(c(COV,"grp"), collapse="+")))
  fitw <- multinom(f, data=d, weights=w, trace=FALSE, maxit=500)
  beta <- matrix(coef(fitw), nrow=2, byrow=FALSE)
  colnames(beta) <- colnames(X)
  vc <- mn_cluster_vcov(fitw, X, as.integer(d$y)-1, w, d$communityID)
  sew <- matrix(sqrt(diag(vc$V)), nrow=2, byrow=FALSE)
  colnames(sew) <- colnames(X)
  cat("clusters used:", vc$G, "\n")
  t_w <- mk_tab(beta, sew, paste0("weighted + cluster-robust (", cut, ")"), nrow(d))
  # unweighted reference
  fitu <- multinom(f, data=d, trace=FALSE, maxit=500)
  bu <- matrix(coef(fitu), nrow=2, byrow=FALSE); colnames(bu) <- colnames(X)
  vu <- mn_cluster_vcov(fitu, X, as.integer(d$y)-1, rep(1,nrow(d)), d$communityID)
  seu <- matrix(sqrt(diag(vu$V)), nrow=2, byrow=FALSE); colnames(seu) <- colnames(X)
  t_u <- mk_tab(bu, seu, paste0("unweighted + cluster-robust (", cut, ")"), nrow(d))
  cat("\n--- weighted multinomial, cluster-robust by community ---\n")
  print(t_w[, c("outcome","group","RRR","lo","hi","p")], row.names=FALSE)
  cat("\n--- unweighted, cluster-robust by community (variance comparison) ---\n")
  print(t_u[, c("outcome","group","RRR","lo","hi","p")], row.names=FALSE)

  # ---- design-based binomials (strata x PSU x weight) ----
  des <- svydesign(ids=~communityID, strata=~stratum, weights=~w13, data=d, nest=TRUE)
  cat("\n--- design-based binomial contrasts (svyglm, strata+PSU+weights) ---\n")
  cat("design df =", degf(des), "\n")
  tb <- list()
  for (oc in c("disabled","dead")) {
    dd <- d; dd$bin <- as.numeric(dd$y==oc); dd <- dd[dd$y %in% c("alive_free", oc), ]
    des2 <- svydesign(ids=~communityID, strata=~stratum, weights=~w13, data=dd, nest=TRUE)
    fb <- as.formula(paste("bin ~", paste(c(COV,"grp"), collapse="+")))
    gw <- svyglm(fb, design=des2, family=quasibinomial())
    co <- summary(gw)$coefficients
    for (g in c("remitted","incident","persistent")) {
      k <- paste0("grp",g)
      est <- co[k,1]; s <- co[k,2]
      tb[[length(tb)+1]] <- data.frame(contrast=paste0(oc," vs alive-free"), group=g,
        RRR=exp(est), lo=exp(est-qt(0.975, degf(des2))*s), hi=exp(est+qt(0.975, degf(des2))*s),
        p=co[k,4], n=nrow(dd), stringsAsFactors=FALSE)
    }
  }
  tb <- do.call(rbind, tb); print(tb, row.names=FALSE)
  res_compare[[cut]] <- list(mn_w=t_w, mn_u=t_u, bin=tb)
  cat("\n")
}

out <- do.call(rbind, lapply(res_compare, function(x) rbind(x$mn_w, x$mn_u)))
write.csv(out, "out/S2_multinom_weighted.csv", row.names=FALSE)
outb <- do.call(rbind, lapply(names(res_compare), function(k){ z <- res_compare[[k]]$bin; z$cutoff <- k; z }))
write.csv(outb, "out/S2_design_binomial.csv", row.names=FALSE)

## ---- sensitivity to the choice of weight ----
cat("--- weight-choice sensitivity (fully adjusted, cut-offs 19) ---\n")
cat("(2018 weights are not usable: participants who died before 2018 have no 2018 weight)\n")
b <- get_base(m, "grp19")
for (wc in c("w11","w13","w13_long")) {
  d <- b[complete.cases(b[, c("y", COV, "grp")]) & !is.na(b[[wc]]) & b[[wc]]>0, ]
  if (length(unique(d$y)) < 3) { cat(wc, ": dropped, outcome has <3 levels\n"); next }
  X <- build_X(d, COV); w <- d[[wc]]/mean(d[[wc]], na.rm=TRUE)
  f <- as.formula(paste("y ~", paste(c(COV,"grp"), collapse="+")))
  fit <- multinom(f, data=d, weights=w, trace=FALSE, maxit=500)
  beta <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(beta) <- colnames(X)
  vc <- mn_cluster_vcov(fit, X, as.integer(d$y)-1, w, d$communityID)
  se <- matrix(sqrt(diag(vc$V)), nrow=2, byrow=FALSE); colnames(se) <- colnames(X)
  cat(sprintf("%-8s n=%d | disabled: remi %.2f inci %.2f pers %.2f | dead: remi %.2f inci %.2f pers %.2f\n",
    wc, nrow(d),
    exp(beta[1,"grpremitted"]), exp(beta[1,"grpincident"]), exp(beta[1,"grppersistent"]),
    exp(beta[2,"grpremitted"]), exp(beta[2,"grpincident"]), exp(beta[2,"grppersistent"])))
}
sink()
cat("\nS2 done -> out/S2_*.csv, out/S2_weighted_report.txt\n")
