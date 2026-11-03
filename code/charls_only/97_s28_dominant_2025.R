# 97_s28_dominant_2025.R -- Supplementary Table S28
# The 2025 age-specific thresholds applied to the AWGS dominant-hand reading
# rather than to the maximum of four trials, in the analytic cohort of 6,052.
# Reviewer M8 asked for this cell; the manuscript previously said it was not reported.
suppressMessages(library(haven))
setwd("E:/workbuddy/JAMDA/analysis")

num <- function(x){ x <- zap_labels(x); if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x) }
conv11 <- function(x) paste0(substr(x, 1, 9), "0", substr(x, 10, 11))
mx <- function(a, b) ifelse(is.na(a), b, ifelse(is.na(b), a, pmax(a, b)))

## ---- trial-level grip and dominant hand -------------------------------------
trials <- function(path, conv = FALSE){
  d <- read_dta(path, col_select = c(ID, qc002, qc003, qc004, qc005, qc006))
  M <- sapply(d[c("qc003","qc004","qc005","qc006")], num)
  M[M %in% c(993, 999)] <- NA
  M[M < 1 | M > 80] <- NA
  L <- mx(M[,1], M[,3]); R <- mx(M[,2], M[,4]); B <- mx(L, R)
  dom <- rep(NA_real_, nrow(M))
  for (k in 1:3){
    w <- which(!is.na(d$qc002) & d$qc002 == k)
    dom[w] <- switch(k, R[w], L[w], B[w])
  }
  id <- as.character(d$ID); if (conv) id <- conv11(id)
  data.frame(ID = id, dom = dom, stringsAsFactors = FALSE)
}
d11 <- trials("_x/2011/biomarkers.dta", TRUE); names(d11)[2] <- "dom11"
d13 <- trials("_x/2013/Biomarker.dta");         names(d13)[2] <- "dom13"

m <- readRDS("out/jamda_full.rds")
m <- merge(m, d11, by = "ID", all.x = TRUE)
m <- merge(m, d13, by = "ID", all.x = TRUE)
coh <- subset(m, a13 %in% 0 & !is.na(grp19) & !is.na(age13) & !is.na(male))
cat(sprintf("analytic cohort n = %d (gate 6,052); dominant hand defined 2011 %d, 2013 %d\n",
            nrow(coh), sum(!is.na(coh$dom11)), sum(!is.na(coh$dom13))))

## ---- cut-offs ---------------------------------------------------------------
low <- function(x, male, age, yr){
  if (yr == 19) ifelse(male == 1, x < 28, ifelse(male == 0, x < 18, NA))
  else          ifelse(male == 1, x < ifelse(age < 65, 34, 28),
                       ifelse(male == 0, x < ifelse(age < 65, 20, 18), NA))
}
coh$m4_19_11  <- low(coh$g11,  coh$male, coh$age13, 19)
coh$m4_19_13  <- low(coh$g13,  coh$male, coh$age13, 19)
coh$m4_25_11  <- low(coh$g11,  coh$male, coh$age11, 25)
coh$m4_25_13  <- low(coh$g13,  coh$male, coh$age13, 25)
coh$dom_19_11 <- low(coh$dom11, coh$male, coh$age13, 19)
coh$dom_19_13 <- low(coh$dom13, coh$male, coh$age13, 19)
coh$dom_25_11 <- low(coh$dom11, coh$male, coh$age11, 25)
coh$dom_25_13 <- low(coh$dom13, coh$male, coh$age13, 25)

cat(sprintf("check against the submitted file: 2011 low %d (file %d), 2013 low %d (file %d)\n",
            sum(coh$m4_19_11 == 1, na.rm = TRUE), sum(coh$low11_19 == 1, na.rm = TRUE),
            sum(coh$m4_19_13 == 1, na.rm = TRUE), sum(coh$low13_19 == 1, na.rm = TRUE)))

## ---- paired contrasts -------------------------------------------------------
summ <- function(d, tag){
  a <- d[[paste0(tag, "_11")]]; b <- d[[paste0(tag, "_13")]]
  p <- !is.na(a) & !is.na(b)
  k11 <- sum(a == 1, na.rm = TRUE); both <- sum(a[p] == 1 & b[p] == 1)
  lo <- a == 1 & !is.na(b)
  data.frame(n = sum(!is.na(a)), prev = 100 * k11 / sum(!is.na(a)),
             rev = 100 * sum(b[lo] == 0, na.rm = TRUE) / sum(lo))
}
boot <- function(d, a11, a13, b11, b13, B = 2000){
  comm <- as.character(d$communityID); u <- unique(comm)
  idx <- split(seq_len(nrow(d)), comm)
  vals <- numeric(B)
  for (i in seq_len(B)){
    pick <- sample(u, length(u), replace = TRUE)
    dd <- d[unlist(idx[pick], use.names = FALSE), ]
    pa <- 100 * mean(dd[[a11]] == 1, na.rm = TRUE); pb <- 100 * mean(dd[[b11]] == 1, na.rm = TRUE)
    ra <- 100 * mean(dd[[a13]][dd[[a11]] == 1] == 0, na.rm = TRUE)
    rb <- 100 * mean(dd[[b13]][dd[[b11]] == 1] == 0, na.rm = TRUE)
    vals[i] <- if (i %% 2 == 0) pb - pa else rb - ra
  }
  list(prev = quantile(vals[seq(2, B, 2)], c(.025, .975)),
       rev  = quantile(vals[seq(1, B, 2)], c(.025, .975)))
}

rows <- list()
for (st in c("whole", "50-64", "65+")){
  d0 <- switch(st, "whole" = coh,
               "50-64" = coh[!is.na(coh$age13) & coh$age13 >= 50 & coh$age13 <= 64, ],
               "65+"   = coh[!is.na(coh$age13) & coh$age13 >= 65, ])
  for (thr in c(25, 19)){
    a11 <- sprintf("m4_%d_11", thr); a13 <- sprintf("m4_%d_13", thr)
    b11 <- sprintf("dom_%d_11", thr); b13 <- sprintf("dom_%d_13", thr)
    p <- d0[!is.na(d0[[a11]]) & !is.na(d0[[b11]]), ]
    r <- p[!is.na(p[[a13]]) & !is.na(p[[b13]]), ]
    ci <- boot(r, a11, a13, b11, b13)
    pa <- 100 * mean(p[[a11]] == 1); pb <- 100 * mean(p[[b11]] == 1)
    ra <- 100 * mean(r[[a13]][r[[a11]] == 1] == 0); rb <- 100 * mean(r[[b13]][r[[b11]] == 1] == 0)
    rows[[length(rows) + 1]] <- data.frame(
      stratum = st, thresholds = thr, n_prev = nrow(p), n_rev = nrow(r),
      prev_max4 = pa, prev_dom = pb, prev_diff = pb - pa, prev_ratio = pb / pa,
      prev_lo = ci$prev[1], prev_hi = ci$prev[2],
      rev_max4 = ra, rev_dom = rb, rev_diff = rb - ra, rev_ratio = rb / ra,
      rev_lo = ci$rev[1], rev_hi = ci$rev[2],
      both_rules = sum(p[[a11]] == 1 & p[[b11]] == 1),
      only_max4 = sum(p[[a11]] == 1 & p[[b11]] == 0),
      only_dom = sum(p[[a11]] == 0 & p[[b11]] == 1))
  }
}
res <- do.call(rbind, rows)
write.csv(res, "out/S28_dominant_2025.csv", row.names = FALSE)
options(width = 200)
print(res, digits = 2, row.names = FALSE)

## ---- kappa between the two reading rules -----------------------------------
kap <- function(a, b){ ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  po <- mean(a == b); pe <- mean(a) * mean(b) + (1 - mean(a)) * (1 - mean(b)); (po - pe) / (1 - pe) }
cat(sprintf("\nagreement between rules (2011): 2019 kappa %.3f, 2025 kappa %.3f\n",
            kap(coh$m4_19_11, coh$dom_19_11), kap(coh$m4_25_11, coh$dom_25_11)))
