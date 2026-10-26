# 20d_sexint.R -- cluster-robust Wald test for sex x exposure interaction (M4 spec, m-test)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")

mn_cov_cl <- function(beta, X, y, w, cluster) {
  K <- 2; p <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P / rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y == k))
  Wm <- list(P[, 2] * (1 - P[, 2]), P[, 3] * (1 - P[, 3]), -P[, 2] * P[, 3])
  bl <- matrix(c(1, 3, 3, 2), 2, 2, byrow = TRUE)
  A <- matrix(0, K * p, K * p)
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[bl[k, l]]]
    A[((k - 1) * p + 1):(k * p), ((l - 1) * p + 1):(l * p)] <- crossprod(X, X * ww)
  }
  S <- matrix(0, N, K * p)
  for (k in 1:K) S[, ((k - 1) * p + 1):(k * p)] <- (w * (Y[, k] - P[, k + 1])) * X
  Sc <- rowsum(S, group = cluster, reorder = FALSE); G <- nrow(Sc)
  Ai <- solve(A)
  (G / (G - 1)) * ((N - 1) / (N - K * p)) * Ai %*% crossprod(Sc) %*% Ai
}

sink("out/20d_sexint.txt", split = TRUE)
for (cut in c("19", "25")) {
  b0 <- get_base(m, paste0("grp", cut))
  d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
    mutate(edu_f = factor(edu), wsc = w13 / mean(w13))
  imp <- mice(d0[, c("age10", "male", "bmi", "nchronic", "cesd", "smoke", "edu_f", "grp", "y")],
              m = 20, seed = 20261023, printFlag = FALSE)
  bs <- list(); Vs <- list()
  for (i in 1:20) {
    di <- complete(imp, i)
    di$edu <- as.numeric(as.character(di$edu_f))
    for (kk in 1:4) di[[paste0("edu", kk)]] <- as.numeric(di$edu == kk)
    di$y <- factor(di$y, levels = c("alive_free", "disabled", "dead"))
    di$grp <- relevel(factor(di$grp, levels = GRP), ref = "never")
    di$communityID <- d0$communityID; di$wsc <- d0$wsc
    X <- model.matrix(~ age10 + male + bmi + nchronic + cesd + smoke + edu1 + edu2 + edu3 + edu4 +
                        grp + male:grp, data = di)
    fit <- multinom(y ~ age10 + male + bmi + nchronic + cesd + smoke + edu1 + edu2 + edu3 + edu4 +
                      grp + male:grp, data = di, weights = di$wsc, trace = FALSE, maxit = 500)
    bb <- matrix(coef(fit), 2, byrow = FALSE); colnames(bb) <- colnames(X)
    V <- mn_cov_cl(bb, X, as.integer(di$y) - 1, di$wsc, di$communityID)
    bs[[i]] <- as.numeric(t(bb)); Vs[[i]] <- V
  }
  mm <- 20; b <- Reduce(`+`, bs) / mm
  W <- Reduce(`+`, Vs) / mm
  Bm <- Reduce(`+`, lapply(bs, function(z) tcrossprod(z - b))) / (mm - 1)
  Tm <- W + (1 + 1 / mm) * Bm
  cn <- colnames(X); p <- length(cn)
  idx <- which(grepl("male:grp", cn))
  ind <- c(idx, p + idx)
  bb <- b[ind]; VV <- Tm[ind, ind]
  st <- as.numeric(t(bb) %*% solve(VV) %*% bb)
  cat(sprintf("\ncut-offs %s: cluster-robust Wald test on %d interaction coefficients\n", cut, length(ind)))
  cat(sprintf("  chi-square = %.2f, df = %d, P = %.4f\n", st, length(ind),
              pchisq(st, length(ind), lower.tail = FALSE)))
  # also the disability-equation-only test
  st1 <- as.numeric(t(b[idx]) %*% solve(Tm[idx, idx]) %*% b[idx])
  cat(sprintf("  disability equation only: chi-square = %.2f, df = %d, P = %.4f\n", st1, length(idx),
              pchisq(st1, length(idx), lower.tail = FALSE)))
  for (k in 1:2) for (g in c("remitted", "incident", "persistent")) {
    j <- which(cn == paste0("grp", g))
    cat(sprintf("    %-8s %-10s RRR %.2f  (interaction coefficient %s)\n",
                c("disabled", "dead")[k], g, exp(b[(k - 1) * p + j]),
                sprintf("%+.3f", b[(k - 1) * p + which(cn == paste0("male:grp", g))])))
  }
}
sink()
