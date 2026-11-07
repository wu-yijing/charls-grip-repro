# jfa_20_cutoff_table.R -- Supplementary Table S41
#
# Two-year reproduction of a low-grip classification under three published
# cut-off sets, in the 6,052-participant landmark sample of the manuscript.
#
# Cut-offs (kg):
#   AWGS 2019   men < 28, women < 18   (primary classification of the manuscript)
#   EWGSOP2     men < 27, women < 16
#   FNIH        men < 26, women < 16   (absolute cut-points; the BMI-adjusted
#                                       variant is not used in this paper)
#
# Everything else is unchanged from the primary analysis: readings are the
# maximum of four trials, both waves are classified with the same cut-off,
# quantities are unweighted, and intervals are percentile cluster-bootstrap
# intervals from 2,000 resamples of communities (seed 20261104), the scheme
# used in Table 1, Table S7, S27 and S31.
#
# Run from the repository root. Input: out/jamda_full.rds (built by 01_build.R).
# Output: out/S41_three_cutoffs.csv and a printed table. The AWGS row must
# reproduce Table 1 (647 low in 2011; 40.5%, 35.6-45.6) -- the script checks it.

suppressMessages(library(dplyr))

d <- readRDS("out/jamda_full.rds")
b <- d %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
stopifnot(nrow(b) == 6052)

g11 <- b$g11; g13 <- b$g13; male <- b$male; cid <- b$communityID
n <- length(g11)

kap <- function(a, bb, cc, dd) {
  N <- a + bb + cc + dd
  po <- (a + dd) / N
  pe <- ((a + bb) * (a + cc) + (cc + dd) * (bb + dd)) / N^2
  (po - pe) / (1 - pe)
}
cnt <- function(l1, l3) {
  t <- table(factor(l1, 0:1), factor(l3, 0:1))
  c(a = t[2, 2], b = t[2, 1], c = t[1, 2], d = t[1, 1])
}

sets <- list(AWGS = c(28, 18), EWGSOP2 = c(27, 16), FNIH = c(26, 16))

byc <- split(seq_len(n), cid); ids <- names(byc)
set.seed(20261104); B <- 2000
boot <- lapply(sets, function(s) matrix(NA_real_, B, 4, dimnames = list(NULL, c("rep", "kappa", "po", "bias"))))

for (s in seq_along(sets)) {
  cut <- ifelse(male == 1, sets[[s]][1], sets[[s]][2])
  l1 <- as.integer(g11 < cut); l3 <- as.integer(g13 < cut)
  for (i in 1:B) {
    pick <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(byc[pick], use.names = FALSE)
    v <- cnt(l1[idx], l3[idx])
    boot[[s]][i, ] <- c(100 * v["a"] / (v["a"] + v["b"]),
                        kap(v[1], v[2], v[3], v[4]),
                        100 * (v["a"] + v["d"]) / sum(v),
                        v[2] / sum(v))
  }
}

out <- list()
for (s in seq_along(sets)) {
  cut <- ifelse(male == 1, sets[[s]][1], sets[[s]][2])
  v <- cnt(as.integer(g11 < cut), as.integer(g13 < cut))
  ci <- apply(boot[[s]], 2, quantile, c(.025, .975), na.rm = TRUE)
  out[[s]] <- data.frame(
    set = names(sets)[s], men = sets[[s]][1], women = sets[[s]][2],
    low11 = sum(v[1:2]), low11p = 100 * sum(v[1:2]) / sum(v),
    low13 = sum(v[c(1, 3)]), low13p = 100 * (v[1] + v[3]) / sum(v),
    rep = 100 * v[1] / (v[1] + v[2]), rep_lo = ci[1, "rep"], rep_hi = ci[2, "rep"],
    kappa = kap(v[1], v[2], v[3], v[4]), k_lo = ci[1, "kappa"], k_hi = ci[2, "kappa"],
    po = 100 * (v[1] + v[4]) / sum(v), pabak = 2 * (v[1] + v[4]) / sum(v) - 1,
    pa = 100 * 2 * v[1] / (2 * v[1] + v[2] + v[3]),
    na = 100 * 2 * v[4] / (2 * v[4] + v[2] + v[3]),
    previdx = abs(v[1] - v[4]) / sum(v),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
print(res, digits = 3)
write.csv(res, "out/S41_three_cutoffs.csv", row.names = FALSE)

# gate: the AWGS row must reproduce Table 1
a <- res[res$set == "AWGS", ]
stopifnot(a$low11 == 647,
          abs(a$rep - 40.5) < 0.05,
          abs(a$rep_lo - 35.6) < 0.2, abs(a$rep_hi - 45.6) < 0.2,
          abs(a$kappa - 0.304) < 0.005)
cat("\nGate passed: AWGS row reproduces Table 1 (647 low in 2011; 40.5%, 35.6-45.6; kappa 0.30).\n")
cat(sprintf("Reproduction across cut-offs: %.1f%% to %.1f%%\n",
            min(res$rep), max(res$rep)))
