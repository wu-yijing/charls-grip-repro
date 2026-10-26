# 20c_table3ci_ipcw.R -- m21 (95% CI for standardized risks) and m9 (IPCW for attrition)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")
fs <- readRDS("out/20b_fit_store.rds")
sink("out/20c_table3ci_ipcw.txt", split = TRUE)

## ---------- m21: standardized risks with 95% CI ----------
cat("===== m21. standardized five-year risks, 95% CI (pooled model, delta method) =====\n")
for (cut in c("19", "25")) {
  r <- fs[[cut]]$rF
  B <- r$B; V <- r$V; cols <- r$cols; p <- r$p
  b0 <- get_base(m, paste0("grp", cut))
  d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
    mutate(wsc = w13 / mean(w13))
  prof <- c(const = 1, age10 = 0, male = 0.5,
            bmi = weighted.mean(d0$bmi, d0$wsc, na.rm = TRUE),
            nchronic = weighted.mean(d0$nchronic, d0$wsc, na.rm = TRUE),
            cesd = weighted.mean(d0$cesd, d0$wsc, na.rm = TRUE),
            smoke = weighted.mean(d0$smoke, d0$wsc, na.rm = TRUE),
            edu1 = 0, edu2 = 0, edu3 = 0, edu4 = 0)
  out <- list()
  for (age in c(60, 65, 70, 75, 80)) for (g in c("never", "remitted", "incident", "persistent")) {
    x <- setNames(rep(0, length(cols)), cols)
    for (v in names(prof)) if (v %in% cols) x[v] <- prof[v]
    x["const"] <- 1; x["male"] <- 0.5; x["age10"] <- (age - 65) / 10
    if (g != "never") x[paste0("grp", g)] <- 1
    e <- as.numeric(B %*% x)
    den <- 1 + exp(e[1]) + exp(e[2])
    p1 <- exp(e[1]) / den; p2 <- exp(e[2]) / den
    g1 <- c(p1 * (1 - p1) * x, -p1 * p2 * x)
    g2 <- c(-p1 * p2 * x, p2 * (1 - p2) * x)
    s1 <- sqrt(as.numeric(t(g1) %*% V %*% g1))
    s2 <- sqrt(as.numeric(t(g2) %*% V %*% g2))
    z <- qt(.975, r$df)
    lg1 <- log(p1 / (1 - p1)); lg2 <- log(p2 / (1 - p2))
    c1 <- plogis(c(lg1 - z * s1 / (p1 * (1 - p1)), lg1 + z * s1 / (p1 * (1 - p1))))
    c2 <- plogis(c(lg2 - z * s2 / (p2 * (1 - p2)), lg2 + z * s2 / (p2 * (1 - p2))))
    out[[length(out) + 1]] <- data.frame(cut = cut, age = age, group = g,
      disabled = 100 * p1, d_lo = 100 * c1[1], d_hi = 100 * c1[2],
      dead = 100 * p2, m_lo = 100 * c2[1], m_hi = 100 * c2[2])
  }
  tab <- do.call(rbind, out)
  tab[, 4:9] <- round(tab[, 4:9], 1)
  print(tab, row.names = FALSE)
  if (cut == "19") write.csv(tab, "out/S6_absolute_risk_ci.csv", row.names = FALSE)
}


mn_cl_vcov <- function(beta, X, y, w, cluster) {
  K <- 2; q <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P / rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y == k))
  Wm <- list(P[, 2] * (1 - P[, 2]), P[, 3] * (1 - P[, 3]), -P[, 2] * P[, 3])
  bl <- matrix(c(1, 3, 3, 2), 2, 2, byrow = TRUE)
  A <- matrix(0, K * q, K * q)
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[bl[k, l]]]
    A[((k - 1) * q + 1):(k * q), ((l - 1) * q + 1):(l * q)] <- crossprod(X, X * ww)
  }
  S <- matrix(0, N, K * q)
  for (k in 1:K) S[, ((k - 1) * q + 1):(k * q)] <- (w * (Y[, k] - P[, k + 1])) * X
  Sc <- rowsum(S, group = cluster, reorder = FALSE); G <- nrow(Sc)
  Ai <- solve(A)
  (G / (G - 1)) * ((N - 1) / (N - K * q)) * Ai %*% crossprod(Sc) %*% Ai
}

## ---------- m9: IPCW sensitivity for differential attrition ----------
cat("\n===== m9. inverse-probability-of-observation weights (complete case) =====\n")
for (cut in c("19", "25")) {
  b0 <- get_base(m, paste0("grp", cut))
  d <- b0 %>% filter(!is.na(w13), w13 > 0) %>% mutate(wsc = w13 / mean(w13))
  d <- d[complete.cases(d[, c(COV, "grp")]), ]
  d$R <- as.numeric(!is.na(d$oc))
  cat(sprintf("\ncut-offs %s: complete cases %d, observed outcome %d\n", cut, nrow(d), sum(d$R)))
  dd <- d[d$R == 1, ]
  dd$y <- factor(dd$y, levels = c("alive_free", "disabled", "dead"))
  obs <- glm(R ~ grp + age10 + male + bmi + nchronic + cesd + smoke,
             data = d, family = binomial())
  d$pobs <- as.numeric(predict(obs, newdata = d, type = "response"))
  cat("  pobs summary:", paste(round(summary(d$pobs[dd$R == 1]), 3), collapse = " "), "\n")
  d$pobs <- predict(obs, newdata = d, type = "response")
  q <- quantile(d$pobs[d$R == 1], c(.01, .99))
  d$ipw <- 1 / pmin(pmax(d$pobs, q[1]), q[2])
  d$ipw <- d$ipw / mean(d$ipw)
  dd2 <- d[d$R == 1, ]
  dd2$y <- factor(dd2$y, levels = c("alive_free", "disabled", "dead"))
  cat("  ipw range:", paste(round(range(dd2$ipw), 3), collapse = "-"), "\n")
  fit <- multinom(y ~ age10 + male + bmi + nchronic + cesd + smoke + edu1 + edu2 + edu3 + edu4 + grp,
                  data = dd2, weights = dd2$wsc * dd2$ipw, trace = FALSE, maxit = 500)
  nm <- colnames(coef(fit)); cf <- matrix(as.numeric(coef(fit)), 2, byrow = FALSE); colnames(cf) <- nm
  Xm <- model.matrix(~ age10 + male + bmi + nchronic + cesd + smoke + edu1 + edu2 + edu3 + edu4 + grp,
                     data = dd2)
  V <- mn_cl_vcov(cf, Xm, as.integer(dd2$y) - 1, dd2$wsc * dd2$ipw, dd2$communityID)
  p <- ncol(cf)
  for (k in 1:2) for (g in c("remitted", "incident", "persistent")) {
    j <- which(nm == paste0("grp", g))
    bb <- cf[k, j]; s <- sqrt(V[(k - 1) * p + j, (k - 1) * p + j])
    cat(sprintf("  %-8s %-10s RRR %.2f (%.2f-%.2f)\n", c("disabled", "dead")[k], g,
                exp(bb), exp(bb - 1.96 * s), exp(bb + 1.96 * s)))
  }
  cat(sprintf("  observation model: median P(observed) %.3f, range %.3f-%.3f; weight truncation at 1st/99th pct\n",
              median(dd$pobs), min(dd$pobs), max(dd$pobs)))
}
sink()
cat("\ndone -> out/20c_table3ci_ipcw.txt, out/S6_absolute_risk_ci.csv\n")
