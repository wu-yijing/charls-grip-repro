# NOT RUNNABLE ON THE SYNTHETIC COHORT -- see charls_only/README.md.
# Reads the 2011 and 2013 CHARLS biomarker files (chair-stand item qh003) to reconcile the
# three published transition analyses; produces Supplementary Table S15. The panel values are
# also archived in out/s15_reconciliation.csv.
.libPaths(c("C:/Users/Administrator/AppData/Local/R/win-library/4.5", .libPaths()))
suppressMessages(library(haven)); suppressMessages(library(dplyr))


m <- readRDS("out/jamda_full.rds") %>% mutate(ID = as.character(ID))
fix11 <- function(x) paste0(substr(x, 1, 10), "0", substr(x, 11, 11))

a <- read_dta("data/charls/2011/biomarkers.dta") %>% mutate(ID = fix11(as.character(ID))) %>%
  transmute(ID, cs1 = suppressWarnings(as.numeric(qh003)))
b <- read_dta("data/charls/2013/Biomarker.dta") %>% mutate(ID = as.character(ID)) %>%
  transmute(ID, cs1b = suppressWarnings(as.numeric(qh003)))

d <- m %>% select(ID, g11, g13, age11, low11_19, low13_19) %>%
  left_join(a, by = "ID") %>% left_join(b, by = "ID") %>%
  filter(!is.na(g11), !is.na(g13), age11 >= 60) %>%
  mutate(cs11 = ifelse(!is.na(cs1) & cs1 > 0 & cs1 < 120, cs1, NA),
         cs13 = ifelse(!is.na(cs1b) & cs1b > 0 & cs1b < 120, cs1b, NA),
         slow11 = !is.na(cs11) & cs11 >= 12,
         slow13 = !is.na(cs13) & cs13 >= 12,
         grip11 = low11_19 == 1,
         grip13 = low13_19 == 1,
         poss11 = grip11 | slow11,
         poss13 = grip13 | slow13)

cat("N (>=60, grip both waves) =", nrow(d), "\n")
cat("chair stand available at both waves =", sum(!is.na(d$cs11) & !is.na(d$cs13)), "\n\n")

g <- d %>% filter(grip11)
cat(sprintf("A1 low grip only:        n=%5d  not reproduced=%5d (%.1f%%)\n",
            nrow(g), sum(!g$grip13), 100 * mean(!g$grip13)))

p <- d %>% filter(poss11)
cat(sprintf("A2 low grip or chair:    n=%5d  not reproduced=%5d (%.1f%%)\n",
            nrow(p), sum(!p$poss13), 100 * mean(!p$poss13)))

p2 <- p %>% mutate(type = case_when(grip11 & slow11 ~ "both",
                                    grip11 ~ "grip only",
                                    TRUE ~ "chair only"))
cat("\nB decomposition (endpoint = no longer possible sarcopenia):\n")
for (k in c("both", "chair only", "grip only")) {
  s <- p2 %>% filter(type == k)
  cat(sprintf("  %-11s n=%5d  %5d (%.1f%%)\n", k, nrow(s), sum(!s$poss13), 100 * mean(!s$poss13)))
}

cat("\nC same cases, endpoint = grip normal in 2013:\n")
for (k in c("both", "chair only", "grip only")) {
  s <- p2 %>% filter(type == k)
  cat(sprintf("  %-11s n=%5d  %5d (%.1f%%)\n", k, nrow(s), sum(!s$grip13), 100 * mean(!s$grip13)))
}

prep <- bind_rows(
  tibble(panel = "A", label = "Low grip only (AWGS 2019 cut-offs)",
         n = nrow(g), rev = sum(!g$grip13), pct = round(100 * mean(!g$grip13), 1)),
  tibble(panel = "A", label = "Low grip or chair stand >= 12 s",
         n = nrow(p), rev = sum(!p$poss13), pct = round(100 * mean(!p$poss13), 1)),
  p2 %>% group_by(type) %>%
    summarise(n = n(), rev = sum(!poss13), pct = round(100 * mean(!poss13), 1), .groups = "drop") %>%
    transmute(panel = "B", label = type, n, rev, pct)
)
print(as.data.frame(prep))
write.csv(prep, "out/s15_reconciliation.csv", row.names = FALSE)
cat("\nwritten: analysis/out/s15_reconciliation.csv\n")
