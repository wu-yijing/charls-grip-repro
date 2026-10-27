# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 73_bounding_mi.R -- bounding analysis under the PRIMARY estimator (weighted, 20-fold MI)
#
# Why this file exists. 70_bounding.R bounded the mortality contrast with complete-case
# covariates. Under that estimator the contrast is 1.88; under the primary (MI) estimator it
# is 2.89. A sensitivity analysis must sit on the same estimator as the estimate it bounds,
# otherwise the range cannot be read against Table 2. The base scenario here reproduces
# Table 2 exactly (1.83 and 2.89), which is the check that the estimator is the right one.
#
# Intervals: Wald intervals from the fitted multinomial, pooled by Rubin's rules. They are
# ~25% narrower than the cluster-robust intervals of Table 2; the footnote of Table S18 says
# so, and no conclusion depends on which is used.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")
sink("out/73_bounding_mi.txt", split = TRUE)

b0 <- get_base(m, "grp19")
base <- b0 %>% filter(!is.na(w13), w13 > 0) %>%
        mutate(edu_f = factor(edu), wsc = w13/mean(w13, na.rm = TRUE))
base$grp <- relevel(factor(base$grp19, levels = GRP), ref = "never")
cat("cohort with both readings =", nrow(base),
    "| without an outcome =", sum(is.na(base$oc)), "\n")
cat("with an observed outcome:", sum(!is.na(base$oc)),
    "| deaths =", sum(base$oc == 2, na.rm = TRUE), "\n\n")

fml <- function(env) as.formula(paste("y ~", paste(c(COV, "grp"), collapse = "+")), env = env)

run_one <- function(dd, seed){
  set.seed(seed)
  imp <- mice(dd[, c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y")],
              m = 20, seed = seed, printFlag = FALSE)
  cells <- paste0(rep(c("disabled","dead"), each = 3), "_",
                  rep(c("remitted","incident","persistent"), 2))
  b <- matrix(NA_real_, 20, 6); s <- matrix(NA_real_, 20, 6)
  colnames(b) <- colnames(s) <- cells
  for (i in 1:20) {
    di <- complete(imp, i)
    di$y <- factor(di$y, levels = c("alive_free","disabled","dead"))
    di$grp <- relevel(factor(di$grp, levels = GRP), ref = "never")
    di$edu <- as.numeric(as.character(di$edu_f))
    for (kk in 1:4) di[[paste0("edu",kk)]] <- as.numeric(di$edu == kk)
    fit <- multinom(fml(environment()), data = di, weights = dd$wsc, trace = FALSE, maxit = 500)
    cf <- summary(fit)$coefficients; se <- summary(fit)$standard.errors
    k <- 0
    for (eq in 1:2) for (g in c("remitted","incident","persistent")) {
      k <- k + 1
      jj <- which(colnames(cf) == paste0("grp", g))
      b[i, k] <- cf[eq, jj]; s[i, k] <- se[eq, jj]
    }
  }
  # Rubin's rules on the log scale
  mb <- colMeans(b); W <- colMeans(s^2)
  Bt <- apply(b, 2, var); Tt <- W + (1 + 1/20) * Bt
  lam <- (1 + 1/20) * Bt / Tt; df <- (20 - 1) / pmax(lam, 1e-12)^2
  df[!is.finite(df) | df <= 0] <- 1e5
  data.frame(cell = cells, RRR = exp(mb),
             lo = exp(mb - qt(.975, df) * sqrt(Tt)),
             hi = exp(mb + qt(.975, df) * sqrt(Tt)),
             n = nrow(dd), row.names = NULL)
}

scen <- list("A. observed outcomes only (base)"          = NULL,
             "B. all without an outcome alive/ADL-free"  = 0,
             "C. all without an outcome with disability" = 1,
             "D. all without an outcome dead"            = 2)

res <- list()
for (nm in names(scen)) {
  dd <- base
  if (is.null(scen[[nm]])) dd <- dd[!is.na(dd$oc), ] else dd$oc[is.na(dd$oc)] <- scen[[nm]]
  dd$y <- factor(dd$oc, levels = 0:2, labels = c("alive_free","disabled","dead"))
  cat("running", nm, "| n =", nrow(dd), "| deaths =", sum(dd$y == "dead"), "\n")
  r <- run_one(dd, seed = 20261027 + which(names(scen) == nm))
  r$scenario <- nm; res[[nm]] <- r
  cat(sprintf("   persistent: disabled %.2f (%.2f-%.2f) | death %.2f (%.2f-%.2f)\n",
              r$RRR[3], r$lo[3], r$hi[3], r$RRR[6], r$lo[6], r$hi[6]))
}
out <- do.call(rbind, res)
write.csv(out, "out/73_bounding_mi.csv", row.names = FALSE)
cat("\n--- persistent vs never-low ---\n")
print(out[grepl("persistent", out$cell), c("scenario","cell","RRR","lo","hi","n")],
      digits = 3, row.names = FALSE)
cat("\n--- all cells, base scenario ---\n")
print(out[out$scenario == names(scen)[1], c("cell","RRR","lo","hi")],
      digits = 3, row.names = FALSE)
