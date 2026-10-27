# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 48_lca.R -- 2-class LCA of three readings under the no-change assumption (identified)
m <- readRDS("out/jamda_full.rds")
sink("out/48_lca.txt", split = TRUE)
b <- subset(m, a13 == 0 & !is.na(grp19) & !is.na(age13) & !is.na(male))
cut <- ifelse(b$male == 1, 28, 18)
Y <- cbind(ifelse(is.na(b$g11), NA, as.integer(b$g11 < cut)),
           ifelse(is.na(b$g13), NA, as.integer(b$g13 < cut)),
           ifelse(is.na(b$g15), NA, as.integer(b$g15 < cut)))
n <- nrow(Y); Tt <- 3

em <- function(Y, maxit = 800, tol = 1e-9) {
  n <- nrow(Y); Tt <- ncol(Y)
  pi <- 0.1; se <- 0.8; sp <- 0.95
  for (it in 1:maxit) {
    lg1 <- rep(0, n); lg0 <- rep(0, n)
    for (t in 1:Tt) {
      o <- !is.na(Y[, t])
      lg1[o] <- lg1[o] + ifelse(Y[o, t] == 1, log(se), log(1 - se))
      lg0[o] <- lg0[o] + ifelse(Y[o, t] == 1, log(1 - sp), log(sp))
    }
    lg1 <- lg1 + log(pi); lg0 <- lg0 + log(1 - pi)
    mx <- pmax(lg1, lg0)
    g <- exp(lg1 - mx) / (exp(lg1 - mx) + exp(lg0 - mx))
    num_s <- den_s <- num_f <- den_f <- 0
    for (t in 1:Tt) {
      o <- !is.na(Y[, t])
      num_s <- num_s + sum(g[o] * Y[o, t]); den_s <- den_s + sum(g[o])
      num_f <- num_f + sum((1 - g[o]) * Y[o, t]); den_f <- den_f + sum(1 - g[o])
    }
    pi_new <- mean(g); se_new <- num_s / den_s; sp_new <- 1 - num_f / den_f
    se_new <- min(max(se_new, 1e-4), 1 - 1e-4); sp_new <- min(max(sp_new, 1e-4), 1 - 1e-4)
    pi_new <- min(max(pi_new, 1e-4), 1 - 1e-4)
    if (max(abs(c(pi_new - pi, se_new - se, sp_new - sp))) < tol) { pi <- pi_new; se <- se_new; sp <- sp_new; break }
    pi <- pi_new; se <- se_new; sp <- sp_new
  }
  # observed vs expected pattern probabilities (complete cases only for G^2)
  pat <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
  pr <- apply(pat, 1, function(y) pi * prod(ifelse(y == 1, se, 1 - se)) + (1 - pi) * prod(ifelse(y == 1, 1 - sp, sp)))
  names(pr) <- apply(pat, 1, paste, collapse = "")
  obs <- table(factor(paste0(Y[,1], Y[,2], Y[,3])))
  obs <- obs[names(pr)]
  G2 <- 2 * sum(obs * log(obs / (pr * sum(obs))))
  list(pi = pi, se = se, sp = sp, pr = pr, obs = obs, G2 = G2, df = 4, g = g, it = it)
}

f <- em(Y)
cat(sprintf("iterations %d\n", f$it))
cat(sprintf("true prevalence of low (no-change model) = %.4f\n", f$pi))
cat(sprintf("sensitivity = %.3f | specificity = %.3f\n", f$se, f$sp))
cat(sprintf("G2 = %.2f on %d df (p = %.4f)\n", f$G2, f$df, pchisq(f$G2, f$df, lower.tail = FALSE)))
cat("\nobserved vs expected (complete cases):\n")
print(round(cbind(obs = as.numeric(f$obs) / sum(f$obs), exp = f$pr), 4))

# observed prevalence, PPV, PPA implied
p_obs <- f$pi * f$se + (1 - f$pi) * (1 - f$sp)
ppa <- (f$pi * f$se^2 + (1 - f$pi) * (1 - f$sp)^2) / p_obs
ppv <- f$pi * f$se / p_obs
cat(sprintf("\nimplied prevalence of a low reading %.3f (observed %.3f)\n", p_obs, mean(Y[,1], na.rm = TRUE)))
cat(sprintf("implied PPA %.3f (observed %.3f) | implied PPV %.3f\n", ppa,
            mean(Y[Y[,1] == 1, 2], na.rm = TRUE), ppv))

# cluster bootstrap for se / sp / pi
set.seed(20261024)
cl <- unique(b$communityID); B <- 200
res <- matrix(NA, B, 3, dimnames = list(NULL, c("pi", "se", "sp")))
for (i in 1:B) {
  pick <- sample(cl, length(cl), replace = TRUE)
  idx <- unlist(lapply(pick, function(k) which(b$communityID == k)))
  ff <- try(em(Y[idx, ], maxit = 500), silent = TRUE)
  if (!inherits(ff, "try-error")) res[i, ] <- c(ff$pi, ff$se, ff$sp)
}
cat("\n--- bootstrap 95% CI (200 cluster resamples) ---\n")
for (j in 1:3) cat(sprintf("%-3s %.3f (%.3f-%.3f)\n", colnames(res)[j],
                           median(res[,j], na.rm=TRUE), quantile(res[,j], .025, na.rm=TRUE),
                           quantile(res[,j], .975, na.rm=TRUE)))
saveRDS(list(fit = f, boot = res), "out/S10_lca.rds")

sink()
