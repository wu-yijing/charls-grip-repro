# 99_stratum_bench_2025.R -- the no-change benchmark under the 2025 cut-offs
#
# Run from the repository root, immediately after 96_stratum_bench.R (it is a
# step in run_all.sh). It reads out/jamda_full.rds, which 01_build.R writes from
# whichever raw files are in data/raw_synthetic; nothing here reads CHARLS data,
# so the pipeline exercises it on the synthetic cohort as a structural check.
# A different analysis directory can be named once with CHARLS_ANALYSIS_DIR; no
# path is hard-coded.
#
# 96_stratum_bench.R runs the benchmark within age strata at the AWGS 2019
# cut-offs. The paper's policy claim, however, is about the 2025 age-specific
# thresholds for the 50-64 band. This script runs the same generative model with
# those cut-offs, so that the 63.7% observed in that stratum can be compared
# with what a population whose strength did not change would produce.
#
# The generative model, the calibration and the error term are unchanged from
# 96_stratum_bench.R: a true 2011 value is drawn from the observed 2011
# distribution, a 2013 reading follows from it with within-person error
# estimated in participants who stayed above the cut-off, and the simulated
# reading is rescaled to the observed marginal. Only the cut-offs and the
# classification change, so the 2019 rows below are a check that this script
# and the published one agree.
#
# The seed and the call order of the first three rows are those of
# 96_stratum_bench.R, including its second calibration call, so that a reader
# who runs both scripts gets the same expectations for the 2019 rows and can
# see that the extension changes nothing but the cut-offs.
#
# The row "50-64, 2025" carries the benchmark quoted in the Results: expected
# 57.3% (95% simulation interval 42.9-71.8) against 63.7% observed, with 433 of
# 3,213 readings low at the 2025 thresholds. The classification itself comes
# from 01_build.R (low11_25 / low13_25, cut-offs <34 and <20 kg in the 50-64
# band, <28 and <18 kg above it).
#
# Reads out/jamda_full.rds. Writes out/99_stratum_bench_2025.csv.

suppressMessages(library(dplyr))

if (!file.exists("out/jamda_full.rds")) {
  alt <- Sys.getenv("CHARLS_ANALYSIS_DIR", "")
  if (nzchar(alt) && dir.exists(alt)) setwd(alt)
}
stopifnot(file.exists("out/jamda_full.rds"))

m <- readRDS("out/jamda_full.rds")
base <- m %>% filter(a13 == 0, !is.na(grp19), !is.na(age13), !is.na(male))

## cut is c(men, women); own = re-estimate the moments within the stratum
sim <- function(dat, cut, calib = "g11", own = FALSE, R = 400) {
  rev <- numeric(0)
  for (sx in c(1, 0)) {
    d <- dat[dat$male == sx, ]
    if (nrow(d) < 30) next
    nv  <- base[base$grp19 == "never" & base$male == sx, ]
    sde <- sd(nv$g13 - nv$g11) / sqrt(2)
    if (own) { mu <- mean(d$g11); sdo <- sd(d$g11) }
    else     { mu <- mean(base$g11[base$male == sx]); sdo <- sd(base$g11[base$male == sx]) }
    r  <- max(0.3, 1 - sde^2 / sdo^2)
    cu <- if (sx == 1) cut[1] else cut[2]
    ref <- if (calib == "g11") d$g11 else d$g13
    l1 <- d$g11 < cu
    if (!any(l1)) next
    for (b in seq_len(R)) {
      Tst <- mu + r * (d$g11 - mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2 * (1 - r), 1e-6)))
      o13 <- Tst + rnorm(nrow(d), 0, sde)
      o13 <- (o13 - mean(o13)) / sd(o13) * sd(ref) + mean(ref)
      rev <- c(rev, 1 - mean((o13 < cu)[l1]))
    }
  }
  rev
}

s5064 <- base[base$age13 >= 50 & base$age13 <= 64, ]
rows <- list()

## ---- 2019 rows, seed and call order from 96_stratum_bench.R ----------------
set.seed(20261024)
for (lab in c("all", "<65", ">=65")) {
  dat <- switch(lab, "all" = base, "<65" = base[base$age13 < 65, ], ">=65" = base[base$age13 >= 65, ])
  e11 <- sim(dat, c(28, 18), "g11")
  e13 <- sim(dat, c(28, 18), "g13")          # run for the RNG stream, as in 96
  l1  <- dat$low11_19 == 1
  obs <- 1 - mean(dat$low13_19[l1])
  rows[[length(rows) + 1]] <- data.frame(
    setting = paste0(lab, ", 2019"), n = nrow(dat), low = sum(l1),
    observed = 100 * obs, expected = 100 * mean(e11),
    lo = 100 * quantile(e11, .025), hi = 100 * quantile(e11, .975),
    ratio = obs / mean(e11), own_moments = FALSE, stringsAsFactors = FALSE)
}

## ---- the 50-64 stratum, both threshold sets --------------------------------
for (own in c(FALSE, TRUE)) {
  for (spec in list(list(tag = "50-64, 2019", cut = c(28, 18), l11 = "low11_19", l13 = "low13_19"),
                    list(tag = "50-64, 2025", cut = c(34, 20), l11 = "low11_25", l13 = "low13_25"))) {
    l1  <- s5064[[spec$l11]] == 1
    obs <- 1 - mean(s5064[[spec$l13]][l1])
    e   <- sim(s5064, spec$cut, "g11", own)
    rows[[length(rows) + 1]] <- data.frame(
      setting = if (own) paste0(spec$tag, " (own moments)") else spec$tag,
      n = nrow(s5064), low = sum(l1), observed = 100 * obs, expected = 100 * mean(e),
      lo = 100 * quantile(e, .025), hi = 100 * quantile(e, .975),
      ratio = obs / mean(e), own_moments = own, stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
# Written through a binary connection so the line endings are LF on every
# platform. R's text connection translates "\n" to CRLF on Windows, and the
# archived copy in manifest/s27/ is stored with LF (see .gitattributes); this
# way the script reproduces that file byte for byte wherever it is run.
con <- file("out/99_stratum_bench_2025.csv", open = "wb")
write.csv(res, con, row.names = FALSE)
close(con)
options(width = 200); print(res, digits = 3, row.names = FALSE)

## ---- closing report -------------------------------------------------------
# The published figures come from the CHARLS cohort (n = 6,052). run_all.sh
# runs this step on the synthetic cohort, where the numbers are a structural
# smoke test; the check below fires only on the real cohort, so a pipeline run
# is never counted as a reproduction of the manuscript's values.
targets <- c("all, 2019" = 59.5, "<65, 2019" = 67.8, ">=65, 2019" = 52.1,
             "50-64, 2025" = 57.3)
n_base <- nrow(base)
cat("\n")
if (n_base == 6052) {
  got <- res$expected[match(names(targets), res$setting)]
  for (i in seq_along(targets)) {
    cat(sprintf("check %-12s expected %5.1f%% | published %5.1f%% | %s\n",
                names(targets)[i], got[i], targets[i],
                if (abs(got[i] - targets[i]) < 0.2) "OK" else "MISMATCH"))
  }
  stopifnot(all(abs(got - targets) < 0.2))
} else {
  cat(sprintf("cohort n = %d (not the CHARLS cohort): the rows above are a structural\n", n_base))
  cat("smoke test. The published run has n = 6,052 and gives 59.5 / 67.8 / 52.1 / 57.3;\n")
  cat("the 57.3% is the benchmark quoted in the Results against 63.7% observed.\n")
}
