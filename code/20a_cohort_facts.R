# 20a_cohort_facts.R -- facts needed for M5, M6, m5, m9, m12 of the rev3 review
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(dplyr)})
m <- readRDS("out/jamda_full.rds")
sink("out/20a_cohort_facts.txt", split = TRUE)

cat("===== M5. Composition of the 3,055 excluded for 'ADL disability in 2013' =====\n")
both <- m %>% filter(!is.na(g11), !is.na(g13), !is.na(low11_19), !is.na(low13_19))
cat("participants with a classifiable grip measurement at both waves:", nrow(both), "\n")
cat("  (of", sum(!is.na(m$g11) & !is.na(m$g13)), "with any grip value at both waves)\n")
cat("\n2013 ADL status among these:\n")
print(table(a13 = both$a13, useNA = "always"))
cat("\n2013 ADL status by 2011 grip group:\n")
print(table(grip2011 = ifelse(both$low11_19 == 1, "low", "normal"), a13 = both$a13, useNA = "always"))
cat("\nRow percentages:\n")
print(round(100 * prop.table(table(ifelse(both$low11_19 == 1, "low", "normal"), both$a13, useNA = "ifany"), 1), 1))

cat("\n-- decomposition of the 9,125 -> 6,070 step --\n")
n_both <- nrow(both)
n_dis  <- sum(both$a13 == 1, na.rm = TRUE)
n_miss <- sum(is.na(both$a13))
n_free <- sum(both$a13 == 0, na.rm = TRUE)
cat(sprintf("total %d = ADL disabled %d + ADL missing %d + ADL-free %d\n", n_both, n_dis, n_miss, n_free))
cat(sprintf("2013 ADL disability rate: %.1f%% ; ADL missing: %.1f%%\n", 100*n_dis/n_both, 100*n_miss/n_both))

cat("\n-- by 2011 grip group: disability, missing, free (n and %) --\n")
for (g in c(0, 1)) {
  s <- both[both$low11_19 == g, ]
  cat(sprintf("2011 %s grip (n=%d): disabled %d (%.1f%%), missing %d (%.1f%%), free %d (%.1f%%)\n",
      ifelse(g == 1, "LOW   ", "normal"),
      nrow(s), sum(s$a13 == 1, na.rm = TRUE), 100*mean(s$a13 == 1, na.rm = TRUE),
      sum(is.na(s$a13)), 100*mean(is.na(s$a13)),
      sum(s$a13 == 0, na.rm = TRUE), 100*mean(s$a13 == 0, na.rm = TRUE)))
}

cat("\n===== m5. 9,125 vs 9,116 =====\n")
cat("grip at both waves (any value):", sum(!is.na(m$g11) & !is.na(m$g13)), "\n")
cat("classifiable at both waves (age/sex known):", n_both, "\n")
cat("difference:", sum(!is.na(m$g11) & !is.na(m$g13)) - n_both, "\n")
cat("of these, missing sex:", sum(!is.na(m$g11) & !is.na(m$g13) & is.na(m$male)),
    "| missing age:", sum(!is.na(m$g11) & !is.na(m$g13) & is.na(m$age11)), "\n")

cat("\n===== kappa: full two-wave sample vs analytic sample =====\n")
tt <- trans_table(m[!is.na(m$g11) & !is.na(m$g13), ], "low11_19", "low13_19")
cat("full two-wave sample: n =", tt$n, " kappa =", round(tt$kappa, 4),
    " observed =", sprintf("%.2f%%", 100*tt$po), " expected =", sprintf("%.2f%%", 100*tt$pe), "\n")
print(tt$tab)
base <- get_base(m, "grp19")
tt2 <- trans_table(base, "low11_19", "low13_19")
cat("analytic sample: n =", tt2$n, " kappa =", round(tt2$kappa, 4),
    " observed =", sprintf("%.2f%%", 100*tt2$po), " expected =", sprintf("%.2f%%", 100*tt2$pe), "\n")
print(tt2$tab)
cat("\nreversion among 2011-low, full sample:", sprintf("%.1f%%", 100*sum(tt$tab[2,1])/sum(tt$tab[2,])), "\n")
cat("reversion among 2011-low, analytic sample:", sprintf("%.1f%%", 100*sum(tt2$tab[2,1])/sum(tt2$tab[2,])), "\n")
cat("incidence among 2011-normal, full sample:", sprintf("%.1f%%", 100*sum(tt$tab[1,2])/sum(tt$tab[1,])), "\n")

cat("\n===== M6. 2011 grip strength, change, and within-person variability =====\n")
b <- base
b$chg <- b$g13 - b$g11
cat("-- weighted means by group (analytic sample, n=6,052) --\n")
for (g in GRP) {
  s <- b[b$grp == g, ]
  cat(sprintf("%-10s n=%4d | g11 %5.1f (%.1f) | g13 %5.1f (%.1f) | change %+5.1f (%.1f)\n",
      g, nrow(s), mean(s$g11, na.rm=TRUE), sd(s$g11, na.rm=TRUE),
      mean(s$g13, na.rm=TRUE), sd(s$g13, na.rm=TRUE),
      mean(s$chg, na.rm=TRUE), sd(s$chg, na.rm=TRUE)))
}
cat("\n-- SD of the within-person difference, never-low group --\n")
nv <- b[b$grp == "never", ]
cat(sprintf("never-low: mean change %+.2f kg, SD of change %.2f kg, n=%d\n",
    mean(nv$chg), sd(nv$chg), nrow(nv)))
cat(sprintf("  paired t-test of zero mean change: t = %.2f, P = %.4g\n",
    t.test(nv$chg)$statistic, t.test(nv$chg)$p.value))
cat(sprintf("  ICC-style: between-person SD of g11 %.2f, of g13 %.2f; within-person SD %.2f\n",
    sd(nv$g11), sd(nv$g13), sd(nv$chg)))
cat(sprintf("  implied measurement-error SD (single reading) = %.2f kg\n", sd(nv$chg)/sqrt(2)))
cat(sprintf("  correlation g11 vs g13 (never-low): r = %.3f\n", cor(nv$g11, nv$g13)))

cat("\n-- change in the two single-reading groups --\n")
for (g in c("remitted", "incident")) {
  s <- b[b$grp == g, ]
  cat(sprintf("%-9s mean change %+5.2f kg (SD %.2f), n=%d\n", g, mean(s$chg), sd(s$chg), nrow(s)))
}

cat("\n-- distance below the cut-off among 2011-low participants --\n")
sd(1)
sdw <- sd(nv$chg)
cut11 <- ifelse(b$male == 1, 28, 18)
b$gap11 <- b$g11 - cut11
low11 <- b[b$low11_19 == 1, ]
cat(sprintf("2011-low participants, n=%d; gap (kg relative to cut-off): median %.1f, IQR %.1f to %.1f\n",
    nrow(low11), median(low11$gap11), quantile(low11$gap11, .25), quantile(low11$gap11, .75)))
cat(sprintf("share within 1 SD of the difference (%.2f kg) below the cut-off: %.1f%%\n",
    sdw, 100*mean(low11$gap11 > -sdw)))
cat(sprintf("share within 2 SD (%.2f kg): %.1f%%\n", 2*sdw, 100*mean(low11$gap11 > -2*sdw)))
cat(sprintf("share within 3 kg of the cut-off: %.1f%%\n", 100*mean(low11$gap11 > -3)))
cat("\nremitted participants only (n=", nrow(b[b$grp=="remitted",]), "):\n", sep="")
rm1 <- b[b$grp == "remitted", ]
cat(sprintf("  gap at 2011: median %.1f kg (IQR %.1f to %.1f); within 1 SD: %.1f%%\n",
    median(rm1$gap11), quantile(rm1$gap11,.25), quantile(rm1$gap11,.75), 100*mean(rm1$gap11 > -sdw)))
cat("\npersistent participants only:\n")
pe1 <- b[b$grp == "persistent", ]
cat(sprintf("  gap at 2011: median %.1f kg (IQR %.1f to %.1f); within 1 SD: %.1f%%\n",
    median(pe1$gap11), quantile(pe1$gap11,.25), quantile(pe1$gap11,.75), 100*mean(pe1$gap11 > -sdw)))

cat("\n-- wave effect: never-low group, 2011 vs 2013 means --\n")
cat(sprintf("2011 %.2f (SE %.3f) vs 2013 %.2f (SE %.3f); difference %+.2f kg (95%% CI %+.2f to %+.2f)\n",
    mean(nv$g11), sd(nv$g11)/sqrt(nrow(nv)), mean(nv$g13), sd(nv$g13)/sqrt(nrow(nv)),
    mean(nv$chg), t.test(nv$chg)$conf.int[1], t.test(nv$chg)$conf.int[2]))

cat("\n===== m9. Loss to follow-up =====\n")
b$lost <- is.na(b$oc)
cat("observed outcome by group:\n")
print(round(100 * prop.table(table(b$grp, b$lost), 1), 1))
cat("\nbaseline characteristics, retained vs lost (whole analytic sample):\n")
for (v in c("age13", "male", "g11", "bmi", "nchronic", "cesd", "smoke")) {
  r <- t.test(b[[v]] ~ b$lost)
  cat(sprintf("%-8s retained %.1f | lost %.1f | diff %+.1f | P = %.4g\n", v,
      mean(b[[v]][!b$lost], na.rm=TRUE), mean(b[[v]][b$lost], na.rm=TRUE),
      mean(b[[v]][!b$lost], na.rm=TRUE) - mean(b[[v]][b$lost], na.rm=TRUE), r$p.value))
}
cat("\nlow-grip prevalence by group and follow-up status (%):\n")
print(round(100 * tapply(b$low13_19, list(b$grp, b$lost), mean), 1))

cat("\n===== m12. Sex-specific prevalence of low grip =====\n")
cat("2011 low-grip prevalence by sex (AWGS 2019 cut-offs, two-wave sample):\n")
for (s in c(0, 1)) {
  ss <- both[both$male == s, ]
  cat(sprintf("  %s: %d/%d = %.1f%%\n", ifelse(s==1,"men  ","women"),
      sum(ss$low11_19), nrow(ss), 100*mean(ss$low11_19)))
}
cat("2011 low-grip prevalence by sex (AWGS 2025 cut-offs):\n")
for (s in c(0, 1)) {
  ss <- m[!is.na(m$low11_25) & m$male == s, ]
  cat(sprintf("  %s: %d/%d = %.1f%%\n", ifelse(s==1,"men  ","women"),
      sum(ss$low11_25), nrow(ss), 100*mean(ss$low11_25)))
}
cat("\nsex composition and low-grip prevalence by group (2019 cut-offs):\n")
for (g in GRP) {
  s <- b[b$grp == g, ]
  cat(sprintf("  %-10s men %.1f%% | 2011 low %.1f%% | 2013 low %.1f%%\n", g,
      100*mean(s$male), 100*mean(s$low11_19), 100*mean(s$low13_19)))
}

cat("\n===== Table 2 denominators for the note (m4) =====\n")
obs <- b[!is.na(b$oc), ]
cat("observed 2018 outcome:", nrow(obs), "| with sampling weight:", sum(!is.na(obs$w13) & obs$w13 > 0), "\n")
cat("complete covariates:", sum(complete.cases(obs[, COV])), "\n")

sink()
cat("\ndone -> out/20a_cohort_facts.txt\n")
