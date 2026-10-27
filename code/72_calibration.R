# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 72_calibration.R -- calibration of the two prediction models compared in the incremental analysis
# Reviewer request (M9): discrimination was reported without calibration, and net benefit depends on it.
# Models are those used for Table 3: binary, complete case, unweighted, competing outcome excluded.
# Reported: apparent Hosmer-Lemeshow and deciles (in sample), plus 10-fold cluster cross-validated
# calibration slope, calibration-in-the-large, Brier score and AUC (out of sample).
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); library(pROC)
m <- readRDS("out/jamda_full.rds")
set.seed(20261027)
sink("out/72_calibration.txt", split = TRUE)

b <- get_base(m, "grp19")
b$g13_5 <- b$g13/5
b$grp4  <- b$grp
d0 <- b[complete.cases(b[, c(COV, "g13_5", "grp4", "oc")]), ]
cat("analysis set n =", nrow(d0), "\n")

hl_test <- function(y, p, g = 10){
  q  <- unique(quantile(p, seq(0, 1, length.out = g + 1)))
  gr <- cut(p, breaks = q, include.lowest = TRUE, labels = FALSE)
  ok <- !is.na(gr)
  o  <- tapply(y[ok], gr[ok], sum); n <- tapply(y[ok], gr[ok], length)
  e  <- tapply(p[ok], gr[ok], sum)
  stat <- sum((o - e)^2 / (e * (1 - e/n)))
  list(stat = stat, df = length(o) - 2, p = 1 - pchisq(stat, length(o) - 2),
       obs = o/n, exp = e/n, n = n)
}

cal <- function(y, p){
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  c(slope = unname(coef(glm(y ~ lp, family = binomial()))[2]),
    citl  = unname(coef(glm(y ~ offset(lp), family = binomial()))[1]),
    brier = mean((y - p)^2),
    auc   = as.numeric(auc(roc(y, p, quiet = TRUE))))
}

cv_pred <- function(dd, fml, K = 10){
  cl <- unique(dd$communityID)
  fold_of <- setNames(sample(rep(1:K, length.out = length(cl))), cl)
  ddp <- rep(NA_real_, nrow(dd))
  for (k in 1:K) {
    tr <- dd[fold_of[dd$communityID] != k, ]
    te <- which(fold_of[dd$communityID] == k)
    fit <- glm(fml, data = tr, family = binomial())
    ddp[te] <- predict(fit, newdata = dd[te, ], type = "response")
  }
  ddp
}

for (outcome in c(1, 2)) {
  lab <- c("disability","death")[outcome]
  dd  <- d0[d0$oc %in% c(0, outcome), ]
  y   <- as.numeric(dd$oc == outcome)
  fA  <- as.formula(paste("y ~", paste(c(COV, "g13_5"), collapse = "+")))
  fB  <- as.formula(paste("y ~", paste(c(COV, "g13_5", "grp4"), collapse = "+")))

  cat("\n===========================================================\n")
  cat(lab, ": n =", nrow(dd), " events =", sum(y),
      " observed risk =", sprintf("%.4f", mean(y)), "\n")
  cat("===========================================================\n")

  for (nm in c("A","B2")) {
    fml <- if (nm == "A") fA else fB
    fit <- glm(fml, data = dd, family = binomial())
    pa  <- predict(fit, type = "response")
    h   <- hl_test(y, pa)
    pc  <- cv_pred(dd, fml)
    cc  <- cal(y, pc)
    ca  <- cal(y, pa)
    cat(sprintf("\n model %-2s  apparent AUC %.3f | CV AUC %.3f | CV Brier %.4f\n", nm, ca["auc"], cc["auc"], cc["brier"]))
    cat(sprintf("           CV calibration slope %.3f | CV CITL %+.4f\n", cc["slope"], cc["citl"]))
    cat(sprintf("           apparent HL %.2f (%d df, P = %.3f) | Emax %.4f\n", h$stat, h$df, h$p, max(abs(h$obs - h$exp))))
    if (nm == "A") {
      cat("\n   deciles, model A (predicted vs observed, n):\n")
      for (i in seq_along(h$obs))
        cat(sprintf("     %2d  %.4f  %.4f  %4d\n", i, h$exp[i], h$obs[i], h$n[i]))
    }
  }
}
sink()
