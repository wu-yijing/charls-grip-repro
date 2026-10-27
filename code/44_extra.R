# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 44_extra.R -- N1 inference + N2 grey zone + N3 threshold sweep + N4 de-attenuation
m <- readRDS("out/jamda_full.rds")
sink("out/44_extra.txt", split = TRUE)
b <- subset(m, a13 == 0 & !is.na(grp19) & !is.na(age13) & !is.na(male))
cut <- ifelse(b$male == 1, 28, 18)
Y <- cbind(ifelse(is.na(b$g11), NA, as.integer(b$g11 < cut)),
           ifelse(is.na(b$g13), NA, as.integer(b$g13 < cut)),
           ifelse(is.na(b$g15), NA, as.integer(b$g15 < cut)))
n <- nrow(Y); Tt <- 3

# ---------- log-likelihood as a function of free parameters ----------
logit <- function(p) log(p / (1 - p)); expit <- function(x) 1 / (1 + exp(-x))

fwd_ll <- function(par) {
  pi2 <- expit(par[1]); se <- expit(par[2]); sp <- expit(par[3])
  a21 <- expit(par[4]); a21b <- expit(par[5])
  pi <- c(1 - pi2, pi2)
  A <- list(matrix(c(1 - a21, a21, 1 - a21, a21), 2, 2, byrow = FALSE),
            matrix(c(1 - a21b, a21b, 1 - a21b, a21b), 2, 2, byrow = FALSE))
  # rebuild rows properly: row j -> transition probs
  A[[1]] <- matrix(c(1 - a21, a21, 1 - a21, a21), 2, 2, byrow = TRUE)
  A[[2]] <- matrix(c(1 - a21b, a21b, 1 - a21b, a21b), 2, 2, byrow = TRUE)
  B <- array(1, dim = c(n, 2, Tt))
  for (t in 1:Tt) {
    obs <- !is.na(Y[, t])
    B[obs, 1, t] <- ifelse(Y[obs, t] == 1, 1 - sp, sp)
    B[obs, 2, t] <- ifelse(Y[obs, t] == 1, se, 1 - se)
  }
  a <- matrix(0, n, 2); for (k in 1:2) a[, k] <- pi[k] * B[, k, 1]
  s <- rowSums(a); ll <- sum(log(pmax(s, 1e-300))); a <- a / s
  for (t in 2:Tt) {
    a <- (a %*% A[[t - 1]]) * B[, , t]
    s <- rowSums(a); ll <- ll + sum(log(pmax(s, 1e-300))); a <- a / s
  }
  ll
}

# MLE from a grid search + optim
best <- NULL
set.seed(20261024)
for (r in 1:12) {
  st <- c(logit(runif(1, .02, .08)), logit(runif(1, .6, .85)), logit(runif(1, .93, .98)),
          logit(runif(1, .05, .3)), logit(runif(1, .05, .3)))
  o <- try(optim(st, fwd_ll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (!inherits(o, "try-error") && (is.null(best) || o$value > best$value)) best <- o
}
par <- best$par; ll <- best$value
se <- expit(par[2]); sp <- expit(par[3]); pi2 <- expit(par[1])
a21 <- expit(par[4]); a21b <- expit(par[5])
cat(sprintf("ll = %.2f\n", ll))
cat(sprintf("se = %.3f (%.3f-%.3f) | sp = %.3f (%.3f-%.3f) | true low wave1 = %.4f | recovery 1: %.3f | recovery 2: %.3f\n",
            se, NA, NA, sp, NA, NA, pi2, a21, a21b))

# ---------- numerical Hessian -> SEs ----------
p <- length(par); H <- matrix(0, p, p); h <- 1e-4
for (i in 1:p) for (j in 1:p) {
  ei <- rep(0, p); ej <- rep(0, p); ei[i] <- h; ej[j] <- h
  H[i, j] <- (fwd_ll(par + ei + ej) - fwd_ll(par + ei - ej) - fwd_ll(par - ei + ej) + fwd_ll(par - ei - ej)) / (4 * h^2)
}
V <- try(solve(-H), silent = TRUE)
if (inherits(V, "try-error")) { V <- MASS::ginv(-H) }
seV <- sqrt(diag(V))
# delta-method CI on the probability scale
ci <- function(k) { est <- expit(par[k]); s <- seV[k] * est * (1 - est); c(est, est - 1.96 * s, est + 1.96 * s) }
tab <- rbind(pi_low = ci(1), sensitivity = ci(2), specificity = ci(3), recovery_1 = ci(4), recovery_2 = ci(5))
cat("\n--- estimates (delta-method 95% CI) ---\n"); print(round(tab, 3))

# ---------- decomposition of observed reversion ----------
emit <- function(state, y, se, sp) ifelse(state == 2, ifelse(y == 1, se, 1 - se), ifelse(y == 1, 1 - sp, sp))
paths <- as.matrix(expand.grid(s1 = 1:2, s2 = 1:2, s3 = 1:2))
pp <- numeric(nrow(paths))
A <- list(matrix(c(1 - a21, a21, 1 - a21, a21), 2, 2, byrow = TRUE),
          matrix(c(1 - a21b, a21b, 1 - a21b, a21b), 2, 2, byrow = TRUE))
for (i in 1:nrow(paths)) { q <- paths[i, ]; pp[i] <- c(1 - pi2, pi2)[q[1]] * A[[1]][q[1], q[2]] * A[[2]][q[2], q[3]] }
pat <- as.matrix(expand.grid(y1 = 0:1, y2 = 0:1, y3 = 0:1))
Ppat <- setNames(numeric(nrow(pat)), apply(pat, 1, paste, collapse = ""))
for (i in 1:nrow(pat)) {
  tot <- 0
  for (s in 1:nrow(paths)) tot <- tot + pp[s] * emit(paths[s,1], pat[i,1], se, sp) * emit(paths[s,2], pat[i,2], se, sp) * emit(paths[s,3], pat[i,3], se, sp)
  Ppat[i] <- tot
}
cat("\nmodel-implied patterns:\n"); print(round(Ppat, 4))
P_y1_low <- Ppat["100"] + Ppat["101"] + Ppat["110"] + Ppat["111"]
P_rev <- (Ppat["100"] + Ppat["101"]) / P_y1_low     # y2 = 0 | y1 = 1
cat(sprintf("implied prevalence of a low 2011 reading = %.3f (observed %.3f)\n", P_y1_low, mean(Y[,1])))
cat(sprintf("implied reversion P(y2=0|y1=1) = %.3f (observed %.3f)\n", P_rev,
            mean(Y[Y[,1]==1,2], na.rm = TRUE)))
# of the observed reversions, how many are true recovery vs misclassification?
# true low at t1 and true normal at t2, with y1=1,y2=0
num_true_rec <- 0
for (s in 1:nrow(paths)) {
  q <- paths[s, ]
  if (q[1] == 2 && q[2] == 1) num_true_rec <- num_true_rec + pp[s] * emit(2,1,se,sp) * emit(1,0,se,sp)
}
cat(sprintf("among apparent reversions: true recovery = %.1f%% ; misclassification = %.1f%%\n",
            100 * num_true_rec / P_rev, 100 * (1 - num_true_rec / P_rev)))

# ---------- N2 grey zone ----------
mdc <- 14.2; half <- mdc / 2
d0 <- b$g11 - cut
inzone <- abs(d0) <= half
low11 <- b$low11_19 == 1
cat("\n--- N2 grey zone (cut-off +/- 7.1 kg) ---\n")
cat(sprintf("whole sample in grey zone: %.1f%% (n = %d)\n", 100 * mean(inzone), sum(inzone)))
cat(sprintf("among 2011 low readings: in zone %.1f%% (n = %d of %d)\n", 100 * mean(inzone[low11]), sum(inzone[low11]), sum(low11)))
rev_in <- mean(b$low13_19[low11 & inzone] == 0); rev_out <- mean(b$low13_19[low11 & !inzone] == 0)
cat(sprintf("reversion: inside zone %.1f%% (n=%d) | outside zone %.1f%% (n=%d)\n",
            100 * rev_in, sum(low11 & inzone), 100 * rev_out, sum(low11 & !inzone)))
# confirmed low (low at both waves) captured by retesting only the zone
conf <- low11 & b$low13_19 == 1
cat(sprintf("confirmed low n = %d; of these %.1f%% lie in the grey zone at 2011\n", sum(conf), 100 * mean(inzone[conf])))
cat(sprintf("retesting only the grey zone = %d tests, capturing %.1f%% of all confirmed lows\n",
            sum(inzone), 100 * sum(conf & inzone) / sum(conf)))
cat(sprintf("retesting all low readings = %d tests\n", sum(low11)))

# ---------- N3 threshold sweep ----------
cat("\n--- N3 threshold sweep (both sexes shifted by the same number of kg) ---\n")
sw <- data.frame()
for (delta in -6:6) {
  cut2 <- cut + delta
  l1 <- as.integer(b$g11 < cut2); l2 <- as.integer(b$g13 < cut2)
  tb <- table(factor(l1, 0:1), factor(l2, 0:1))
  po <- sum(diag(tb)) / sum(tb); pe <- sum(rowSums(tb) * colSums(tb)) / sum(tb)^2
  kap <- (po - pe) / (1 - pe)
  ppa <- 2 * tb[2,2] / (2 * tb[2,2] + tb[2,1] + tb[1,2])
  rev <- tb[2,1] / sum(tb[2,])
  sw <- rbind(sw, data.frame(delta = delta, cut_men = 28 + delta, prev = mean(l1), kappa = kap,
                             ppa = ppa, reversion = rev))
}
print(round(sw, 3), row.names = FALSE)
write.csv(sw, "out/S11_threshold_sweep.csv", row.names = FALSE)

# ---------- N4 de-attenuation ----------
cat("\n--- N4 regression-dilution correction ---\n")
for (R in c(0.67, 0.75, 0.80)) {
  cat(sprintf("R = %.2f : disability 1.83 -> %.2f (%.2f-%.2f) | death 2.89 -> %.2f (%.2f-%.2f)\n", R,
              exp(log(1.83)/R), exp(log(1.14)/R), exp(log(2.94)/R),
              exp(log(2.89)/R), exp(log(1.73)/R), exp(log(4.83)/R)))
}

out <- list(tab = tab, Ppat = Ppat, se = se, sp = sp, pi2 = pi2, a21 = a21, a21b = a21b,
            rev_in = rev_in, rev_out = rev_out, sweep = sw, ll = ll)
saveRDS(out, "out/S10_extra.rds")
cat("\nsaved out/S10_extra.rds\n")

sink()
