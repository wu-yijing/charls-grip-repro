source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr)})
m <- readRDS("out/jamda_full.rds")
b0 <- get_base(m, "grp19")
d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
  mutate(nr = low11_19 + low13_19, edu_f = factor(edu), wsc = w13 / mean(w13))
d0$nr_low <- d0$nr
imp <- mice(d0[, c("age10", "male", "bmi", "nchronic", "cesd", "smoke", "edu_f", "grp", "y", "nr_low")],
            m = 20, seed = 20261023, printFlag = FALSE)
bs <- list(); vs <- list()
for (i in 1:20) {
  di <- complete(imp, i)
  di$y <- factor(di$y, levels = c("alive_free", "disabled", "dead"))
  di$grp <- relevel(factor(di$grp, levels = GRP), ref = "never")
  di$edu <- as.numeric(as.character(di$edu_f))
  for (k in 1:4) di[[paste0("edu", k)]] <- as.numeric(di$edu == k)
  f <- as.formula(paste("y ~", paste(c(COV, "grp"), collapse = "+")))
  fit <- multinom(f, data = di, weights = d0$wsc, trace = FALSE, maxit = 500)
  cf <- summary(fit)$coefficients; se <- summary(fit)$standard.errors
  bs[[i]] <- as.numeric(t(cf)); vs[[i]] <- as.numeric(t(se^2))
}
b <- Reduce(`+`, bs) / 20; W <- Reduce(`+`, vs) / 20
Bv <- Reduce(`+`, lapply(bs, function(z) (z - b)^2)) / 19; Tm <- W + (1 + 1/20) * Bv
p <- length(b) / 2
for (k in 1:2) for (g in c("remitted", "incident", "persistent")) {
  j <- which(colnames(summary(fit)$coefficients) == paste0("grp", g)); idx <- (k - 1) * p + j
  cat(sprintf("%-9s %-11s %.3f (%.2f-%.2f)\n", c("disabled", "dead")[k], g,
              exp(b[idx]), exp(b[idx] - qt(.975, 356) * sqrt(Tm[idx])),
              exp(b[idx] + qt(.975, 356) * sqrt(Tm[idx]))))
}
