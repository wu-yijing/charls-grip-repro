# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 74_frailty_mi.R -- frailty index excluding grip, under the PRIMARY estimator (weighted MI)
#
# 71_frailty.R fitted the index on complete cases, where the mortality contrast is 1.88.
# Under the primary estimator it is 2.89, so the "1.88 -> 1.87" comparison was being read on
# the wrong scale. This script repeats the analysis on the MI estimator so that the two
# numbers sit next to Table 2 rather than next to a complete-case side analysis.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(haven); library(dplyr)})
m <- readRDS("out/jamda_full.rds")
sink("out/74_frailty_mi.txt", split = TRUE)

num <- function(x){ x <- zap_labels(x); if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x) }
conv11 <- function(x) paste0(substr(x,1,9), "0", substr(x,10,11))

CHRONIC <- paste0("da007_", 1:14, "_")
h11 <- read_dta("data/raw_synthetic/2011/health_status_and_functioning.dta", col_select = all_of(c("ID", CHRONIC)))
h11 <- data.frame(ID = conv11(as.character(h11$ID)), h11[CHRONIC])
ch <- sapply(h11[CHRONIC], function(s) as.numeric(grepl("^1", as.character(num(s)))))
h11 <- data.frame(ID = h11$ID, ch)

CESD <- c("dc009","dc010","dc011","dc012","dc013","dc014","dc015","dc016","dc017","dc018")
h13 <- read_dta("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", col_select = all_of(c("ID","da002", CESD)))
h13 <- data.frame(ID = as.character(h13$ID), srh = num(h13$da002), h13[CESD])
Cm <- sapply(h13[CESD], num); Cm[, c("dc013","dc016")] <- 5 - Cm[, c("dc013","dc016")]
ced <- ifelse(is.na(Cm), NA, as.numeric(Cm >= 3)); colnames(ced) <- paste0("cesd_", 1:10)
h13 <- data.frame(ID = h13$ID, srh = h13$srh, ced)
h13$srh_def <- ifelse(is.na(h13$srh), NA, as.numeric(h13$srh >= 4))

fi <- m %>% select(ID, bmi) %>% mutate(bmi_def = ifelse(is.na(bmi), NA, as.numeric(bmi < 18.5)))
fi <- fi %>% left_join(h11, by = "ID") %>% left_join(h13, by = "ID")
ITEMS <- c(CHRONIC, paste0("cesd_", 1:10), "srh_def", "bmi_def")
X <- as.matrix(fi[, ITEMS]); nobs <- rowSums(!is.na(X)); ntot <- rowSums(X, na.rm = TRUE)
fi$fi <- ifelse(nobs >= 0.8 * length(ITEMS), ntot / nobs, NA)

b0 <- get_base(m, "grp19")
base <- b0 %>% filter(!is.na(w13), w13 > 0, !is.na(oc)) %>%
        left_join(fi %>% select(ID, fi), by = "ID") %>%
        mutate(edu_f = factor(edu), wsc = w13/mean(w13, na.rm = TRUE))
base$grp <- relevel(factor(base$grp19, levels = GRP), ref = "never")
base$y   <- factor(base$oc, levels = 0:2, labels = c("alive_free","disabled","dead"))
cat("analytic n =", nrow(base), "| FI available =", sum(!is.na(base$fi)),
    "| deaths =", sum(base$oc == 2), "\n")
cat("\nFI by group (mean, SD):\n")
print(base %>% group_by(grp19) %>% summarise(n = n(),
      FI = sprintf("%.3f (%.3f)", mean(fi, na.rm = TRUE), sd(fi, na.rm = TRUE)),
      missing = sum(is.na(fi))))
cat("\nFI tertiles (whole cohort):", round(quantile(base$fi, c(1/3, 2/3), na.rm = TRUE), 3), "\n")
base$fi_t <- cut(base$fi, breaks = c(-Inf, quantile(base$fi, c(1/3, 2/3), na.rm = TRUE), Inf),
                 labels = c("low","middle","high"))

pool_fit <- function(dd, extra, seed, label){
  set.seed(seed)
  vars <- c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y", extra)
  imp <- mice(dd[, vars], m = 20, seed = seed, printFlag = FALSE)
  res <- list()
  for (i in 1:20) {
    di <- complete(imp, i)
    di$y <- factor(di$y, levels = c("alive_free","disabled","dead"))
    di$grp <- relevel(factor(di$grp, levels = GRP), ref = "never")
    di$edu <- as.numeric(as.character(di$edu_f))
    for (kk in 1:4) di[[paste0("edu",kk)]] <- as.numeric(di$edu == kk)
    f <- as.formula(paste("y ~", paste(c(COV, "grp", extra), collapse = "+")), env = environment())
    fit <- multinom(f, data = di, weights = dd$wsc, trace = FALSE, maxit = 500)
    cf <- summary(fit)$coefficients; se <- summary(fit)$standard.errors
    jp <- which(colnames(cf) == "grppersistent")
    out <- data.frame(eq = c("disabled","dead"),
                      b_p = cf[, jp], s_p = se[, jp],
                      b_fi = NA_real_, s_fi = NA_real_)
    jf <- which(colnames(cf) == "fi")
    if (length(jf) == 1) { out$b_fi <- cf[, jf]; out$s_fi <- se[, jf] }
    out$imp <- i; res[[i]] <- out
  }
  R <- do.call(rbind, res)
  agg <- R %>% group_by(eq) %>%
    summarise(b = mean(b_p), s = sqrt(mean(s_p^2) + (1 + 1/20) * var(b_p)),
              bf = mean(b_fi), sf = sqrt(mean(s_fi^2, na.rm = TRUE) + (1 + 1/20) * var(b_fi)),
              .groups = "drop")
  data.frame(model = label, outcome = agg$eq,
             RRR = exp(agg$b), lo = exp(agg$b - 1.96*agg$s), hi = exp(agg$b + 1.96*agg$s),
             fi_per_0.1 = exp(agg$bf * 0.1),
             fi_lo = exp((agg$bf - 1.96*agg$sf) * 0.1),
             fi_hi = exp((agg$bf + 1.96*agg$sf) * 0.1), n = nrow(dd))
}

cat("\n--- models (weighted MI, 20 imputations) ---\n")
r1 <- pool_fit(base, character(0), 20261031, "without the index (primary)")
r2 <- pool_fit(base, "fi",        20261032, "with the index, per 0.1")
r3 <- pool_fit(base, "fi_t",      20261033, "with the index in tertiles")
out <- rbind(r1, r2, r3)
print(out[, c("model","outcome","RRR","lo","hi","fi_per_0.1","fi_lo","fi_hi","n")],
      digits = 3, row.names = FALSE)
write.csv(out, "out/74_frailty_mi.csv", row.names = FALSE)
