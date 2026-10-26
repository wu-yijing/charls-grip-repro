# 20b_design_mi.R -- M3 / M4 / M8 / m11 / m21 of the rev3 review
# (1) design-based multinomial logit with the SAME estimand as the primary model, with MI
# (2) between-group Wald contrasts and an ordinal trend test
# (3) one table of every sensitivity specification (Supplementary Table S5)
# (4) 95% CIs for the standardized risks (Table 3 / new Supplementary Table S6)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")

strat_cluster_var <- function(Z, cluster, stratum) {
  Zt <- rowsum(Z, group = cluster, reorder = FALSE)
  st <- stratum[match(rownames(Zt), as.character(cluster))]
  gm <- colMeans(Zt)
  V <- matrix(0, ncol(Z), ncol(Z))
  for (h in unique(st)) {
    idx <- which(st == h); nh <- length(idx)
    if (nh >= 2) {
      Zh <- sweep(Zt[idx, , drop = FALSE], 2, colMeans(Zt[idx, , drop = FALSE]))
      V <- V + (nh / (nh - 1)) * crossprod(Zh)
    } else {
      d <- Zt[idx, ] - gm
      V <- V + tcrossprod(d)
    }
  }
  V
}

# multinomial score / bread; returns both cluster-robust and design-based covariance
mn_cov <- function(beta, X, y, w, cluster, stratum) {
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
  Ai <- solve(A)
  Sc <- rowsum(S, group = cluster, reorder = FALSE); G <- nrow(Sc)
  Vcl <- (G / (G - 1)) * ((N - 1) / (N - K * p)) * Ai %*% crossprod(Sc) %*% Ai
  Vdb <- Ai %*% strat_cluster_var(S, cluster, stratum) %*% Ai
  list(Vcl = Vcl, Vdb = Vdb)
}

rubin_full <- function(bs, Vs) {
  mm <- length(bs); b <- Reduce(`+`, bs) / mm
  W <- Reduce(`+`, Vs) / mm
  Bm <- Reduce(`+`, lapply(bs, function(z) tcrossprod(z - b))) / (mm - 1)
  T <- W + (1 + 1 / mm) * Bm
  lam <- mean(diag((1 + 1 / mm) * Bm) / diag(T))
  df <- (mm - 1) / max(lam, 1e-12)^2
  list(b = b, V = T, df = max(df, 1))
}

fit_spec <- function(dat, covars, expo = "grp", weights = TRUE, vartype = "cl",
                     outcome = "y", extra = NULL) {
  bs <- list(); Vs <- list(); cn <- NULL
  for (i in seq_along(dat)) {
    di <- dat[[i]]
    if (is.numeric(di[[outcome]]))
      di[[outcome]] <- factor(di[[outcome]], levels = 0:2, labels = c("alive_free", "disabled", "dead"))
    else
      di[[outcome]] <- factor(di[[outcome]], levels = c("alive_free", "disabled", "dead"))
    di$grp <- relevel(factor(di$grp, levels = GRP), ref = "never")
    X <- cbind(1, as.matrix(di[, covars]))
    if (expo == "grp") for (g in c("remitted", "incident", "persistent")) X <- cbind(X, as.numeric(di$grp == g))
    if (expo == "nr") X <- cbind(X, di$nr_low)
    cn <- c("const", covars, if (expo == "grp") paste0("grp", c("remitted", "incident", "persistent")) else "nr_low")
    colnames(X) <- cn
    w <- if (weights) di$wsc else rep(1, nrow(di))
    rhs <- paste(c(covars, if (expo == "grp") "grp" else "nr_low"), collapse = " + ")
    f <- as.formula(paste(outcome, "~", rhs))
    fit <- multinom(f, data = di, weights = w, trace = FALSE, maxit = 500, reltol = 1e-10)
    bb <- matrix(coef(fit), nrow = 2, byrow = FALSE); colnames(bb) <- cn
    vv <- mn_cov(bb, X, as.integer(di[[outcome]]) - 1, w, di$communityID, di$stratum)
    bs[[i]] <- as.numeric(t(bb)); Vs[[i]] <- if (vartype == "db") vv$Vdb else vv$Vcl
  }
  pr <- rubin_full(bs, Vs); p <- length(cn)
  list(B = matrix(pr$b, 2, byrow = TRUE), V = pr$V, df = pr$df, cols = cn, p = p)
}

rows_of <- function(res, label, spec) {
  out <- list()
  for (k in 1:2) for (g in c("remitted", "incident", "persistent")) {
    j <- which(res$cols == paste0("grp", g))
    idx <- (k - 1) * res$p + j
    b <- res$B[k, j]; s <- sqrt(res$V[idx, idx]); df <- res$df
    out[[length(out) + 1]] <- data.frame(
      spec = spec, label = label, outcome = c("disability", "death")[k], group = g,
      est = exp(b), lo = exp(b - qt(.975, df) * s), hi = exp(b + qt(.975, df) * s),
      p = 2 * pt(-abs(b / s), df), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

sink("out/20b_design_mi.txt", split = TRUE)
cat("===== M3/M4/M8: design-based multinomial, Wald contrasts, S5 specifications =====\n")

all_specs <- list(); wald <- list(); trend <- list(); fit_store <- list()

for (cut in c("19", "25")) {
  grpcol <- paste0("grp", cut)
  b0 <- get_base(m, grpcol)
  d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
    mutate(nr_low = if (cut == "19") low11_19 + low13_19 else low11_25 + low13_25,
           edu_f = factor(edu), wsc = w13 / mean(w13))
  cat(sprintf("\n### cut-offs %s | n = %d | clusters = %d | strata = %d\n",
              cut, nrow(d0), length(unique(d0$communityID)), length(unique(d0$stratum))))

  imp <- mice(d0[, c("age10", "male", "bmi", "nchronic", "cesd", "smoke", "edu_f",
                     "grp", "y", "nr_low")],
              m = 20, seed = 20261023, printFlag = FALSE)
  dats <- lapply(1:20, function(i) {
    di <- complete(imp, i)
    di$edu <- as.numeric(as.character(di$edu_f))
    for (kk in 1:4) di[[paste0("edu", kk)]] <- as.numeric(di$edu == kk)
    di$grp <- factor(di$grp, levels = GRP)
    di$y <- di$oc <- as.numeric(di$y) - 1
    di$y4 <- d0$oc4; di$yC <- d0$ocC
    di$communityID <- d0$communityID; di$stratum <- d0$stratum; di$wsc <- d0$wsc
    di
  })
  dats_nw <- lapply(dats, function(di) { di$wsc <- 1; di })
  dats_nobmi <- lapply(dats, function(di) { di$bmi <- 0; di })

  rA  <- fit_spec(dats, c("age10", "male"), "grp", TRUE, "cl")
  rF  <- fit_spec(dats, COV, "grp", TRUE, "cl")
  rDb <- fit_spec(dats, COV, "grp", TRUE, "db")
  rUw <- fit_spec(dats_nw, COV, "grp", FALSE, "cl")
  rT  <- fit_spec(dats, COV, "nr", TRUE, "cl")
  rNB <- fit_spec(dats_nobmi, c("age10", "male", "nchronic", "cesd", "smoke", "edu1", "edu2", "edu3", "edu4"),
                  "grp", TRUE, "cl")
  rD4 <- fit_spec(dats, COV, "grp", TRUE, "cl", outcome = "y4")
  rDC <- fit_spec(dats, COV, "grp", TRUE, "cl", outcome = "yC")

  all_specs[[paste0(cut, "_ageSex")]] <- rows_of(rA, "Age- and sex-adjusted, weighted MI", "A")
  all_specs[[paste0(cut, "_prim")]]   <- rows_of(rF, "Primary: weighted MI, cluster-robust", "B")
  all_specs[[paste0(cut, "_dbmn")]]   <- rows_of(rDb, "Design-based multinomial, weighted MI", "C")
  all_specs[[paste0(cut, "_uw")]]     <- rows_of(rUw, "Unweighted MI, cluster-robust", "D")
  all_specs[[paste0(cut, "_nobmi")]]  <- rows_of(rNB, "Weighted MI, BMI omitted", "E")
  all_specs[[paste0(cut, "_strict")]] <- rows_of(rD4, "Strict disability (inability only), weighted MI", "F")
  all_specs[[paste0(cut, "_comp")]]   <- rows_of(rDC, "Composite outcome incl. IADL, weighted MI", "G")
  fit_store[[cut]] <- list(rF = rF, rDb = rDb, rT = rT)

  cat("\n-- primary, weighted MI, cluster-robust --\n"); print(rows_of(rF, "", "B")[, c("outcome","group","est","lo","hi","p")], row.names = FALSE, digits = 3)
  cat("\n-- design-based multinomial (strata + PSU + weights), weighted MI --\n"); print(rows_of(rDb, "", "C")[, c("outcome","group","est","lo","hi","p")], row.names = FALSE, digits = 3)
  cat("\n-- other specifications --\n")
  for (nm in c("ageSex", "uw", "nobmi", "strict", "comp")) {
    cat("  [", nm, "]\n")
    print(all_specs[[paste0(cut, "_", nm)]][, c("outcome","group","est","lo","hi")], row.names = FALSE, digits = 3)
  }
  cat("\n-- ordinal trend, per additional low reading --\n")
  for (k in 1:2) {
    j <- which(rT$cols == "nr_low"); idx <- (k - 1) * rT$p + j
    bb <- rT$B[k, j]; s <- sqrt(rT$V[idx, idx])
    cat(sprintf("   %-10s RRR %.2f (%.2f-%.2f) P = %.4g\n", c("disability","death")[k],
                exp(bb), exp(bb - qt(.975, rT$df) * s), exp(bb + qt(.975, rT$df) * s),
                2 * pt(-abs(bb / s), rT$df)))
    trend[[paste0(cut, k)]] <- data.frame(cut = cut, outcome = c("disability","death")[k],
      RRR = exp(bb), lo = exp(bb - qt(.975, rT$df) * s), hi = exp(bb + qt(.975, rT$df) * s),
      p = 2 * pt(-abs(bb / s), rT$df))
  }

  # between-group Wald contrasts, within the primary model
  cat("\n-- between-group Wald contrasts (primary model) --\n")
  cmp <- list(c("remitted", "incident"), c("incident", "persistent"), c("remitted", "persistent"))
  for (k in 1:2) for (cc in cmp) {
    j1 <- which(rF$cols == paste0("grp", cc[1])); j2 <- which(rF$cols == paste0("grp", cc[2]))
    i1 <- (k - 1) * rF$p + j1; i2 <- (k - 1) * rF$p + j2
    d <- rF$B[k, j1] - rF$B[k, j2]
    s <- sqrt(rF$V[i1, i1] + rF$V[i2, i2] - 2 * rF$V[i1, i2])
    cat(sprintf("   %-10s %s vs %s: ratio %.2f (%.2f-%.2f), Wald P = %.3f\n",
                c("disability","death")[k], cc[1], cc[2], exp(d),
                exp(d - qt(.975, rF$df) * s), exp(d + qt(.975, rF$df) * s),
                2 * pt(-abs(d / s), rF$df)))
    wald[[paste(cut, k, cc[1], cc[2], sep = "_")]] <- data.frame(cut = cut,
      outcome = c("disability","death")[k], contrast = paste(cc[1], "vs", cc[2]),
      ratio = exp(d), lo = exp(d - qt(.975, rF$df) * s), hi = exp(d + qt(.975, rF$df) * s),
      p = 2 * pt(-abs(d / s), rF$df))
  }

  # sex interaction (LR test on each imputation, then median P)
  cat("\n-- sex x exposure interaction --\n")
  ps6 <- ps3 <- numeric(20)
  for (i in 1:20) {
    di <- dats[[i]]
    di$y <- factor(di$y, levels = 0:2, labels = c("alive_free", "disabled", "dead"))
    f0 <- as.formula(paste("y ~", paste(c(COV, "grp"), collapse = "+")))
    f1 <- as.formula(paste("y ~", paste(c(COV, "grp", "male:grp"), collapse = "+")))
    l0 <- logLik(multinom(f0, data = di, weights = di$wsc, trace = FALSE, maxit = 500))
    l1 <- logLik(multinom(f1, data = di, weights = di$wsc, trace = FALSE, maxit = 500))
    st <- 2 * (as.numeric(l1) - as.numeric(l0))
    ps6[i] <- pchisq(st, df = 6, lower.tail = FALSE)
    ps3[i] <- pchisq(st, df = 3, lower.tail = FALSE)
  }
  cat(sprintf("   median P (df = 6, two equations) = %.3f (range %.3f-%.3f)\n", median(ps6), min(ps6), max(ps6)))
  cat(sprintf("   median P (df = 3, disability equation only) = %.3f\n", median(ps3)))
}
sink()

s5 <- do.call(rbind, all_specs)
write.csv(s5, "out/S5_sensitivity_all.csv", row.names = FALSE)
write.csv(do.call(rbind, wald), "out/S5_wald.csv", row.names = FALSE)
write.csv(do.call(rbind, trend), "out/S5_trend.csv", row.names = FALSE)
saveRDS(fit_store, "out/20b_fit_store.rds")
cat("\nwritten out/S5_sensitivity_all.csv, S5_wald.csv, S5_trend.csv\n")
