# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 40_share_vs_charls.R -- where does the AWGS cut-off sit in CHARLS vs SHARE?
suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
sink("out/40_share_vs_charls.txt", split = TRUE)

d <- m %>%
  filter(!is.na(g11), !is.na(age13), !is.na(male)) %>%
  mutate(age11 = age13 - 2,
         band = cut(age11, breaks = c(49, 54, 59, 64, 69, 74, 79, 84, 120),
                    labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")),
         sex = ifelse(male == 1, "men", "women")) %>%
  filter(age11 >= 50)

cat("CHARLS 2011, age 50+, valid grip: n =", nrow(d), "\n")
cat("men:", sum(d$sex == "men"), " women:", sum(d$sex == "women"), "\n\n")

qs <- c(.05,.10,.20,.30,.40,.50,.60,.70,.80,.90,.95)
for (s in c("men","women")) {
  cat("=== ", s, " ===\n")
  cat(sprintf("%-7s %5s %s\n", "band", "n", paste(sprintf("P%-4d", qs*100), collapse=" ")))
  for (b in levels(d$band)) {
    x <- d$g11[d$sex == s & d$band == b]
    if (length(x) < 30) { cat(sprintf("%-7s %5d  (n too small)\n", b, length(x))); next }
    cat(sprintf("%-7s %5d %s\n", b, length(x),
        paste(sprintf("%-5.0f", quantile(x, qs, na.rm=TRUE, type=7)), collapse=" ")))
  }
  x <- d$g11[d$sex == s]
  cat(sprintf("%-7s %5d %s\n", "pooled", length(x),
      paste(sprintf("%-5.0f", quantile(x, qs, na.rm=TRUE, type=7)), collapse=" ")))
  # where does the AWGS cut-off sit?
  cut_kg <- ifelse(s == "men", 28, 18)
  cat(sprintf("  empirical percentile of the AWGS cut-off (%g kg): %.1f%%\n", cut_kg,
              100 * mean(x < cut_kg)))
  cat(sprintf("  prevalence of low grip: %.1f%%\n\n", 100 * mean(x < cut_kg)))
}

# by band: percentile position of the cut-off
cat("=== percentile position of the AWGS cut-off, by age band ===\n")
for (s in c("men","women")) {
  cut_kg <- ifelse(s == "men", 28, 18)
  for (b in levels(d$band)) {
    x <- d$g11[d$sex == s & d$band == b]
    if (length(x) < 30) next
    cat(sprintf("%-5s %-7s n=%4d  cut=%2d kg -> P%4.1f  (low-grip prevalence %.1f%%)\n",
        s, b, length(x), cut_kg, 100*mean(x < cut_kg), 100*mean(x < cut_kg)))
  }
}

sink()
