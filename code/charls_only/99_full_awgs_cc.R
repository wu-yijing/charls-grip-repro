# 99_full_awgs_cc.R -- chair-stand complete-case check for Table S37 (v1.0.8)
# Run from the analysis directory described in code/charls_only/README.md.
suppressMessages({library(haven); library(dplyr)})
m <- readRDS("out/jamda_full.rds") %>% mutate(ID = as.character(ID))
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
fix11 <- function(x) paste0(substr(x,1,10), "0", substr(x,11,11))
a <- read_dta("_x/2011/biomarkers.dta") %>% mutate(ID = fix11(as.character(ID))) %>%
  transmute(ID, cs1 = suppressWarnings(as.numeric(qh003)))
c13 <- read_dta("_x/2013/Biomarker.dta") %>% mutate(ID = as.character(ID)) %>%
  transmute(ID, cs1b = suppressWarnings(as.numeric(qh003)))
d <- b %>% select(ID, low11_19, low13_19, communityID) %>%
  left_join(a, by="ID") %>% left_join(c13, by="ID") %>%
  mutate(cs11 = ifelse(!is.na(cs1) & cs1 > 0 & cs1 < 120, cs1, NA),
         cs13 = ifelse(!is.na(cs1b) & cs1b > 0 & cs1b < 120, cs1b, NA),
         slow11 = !is.na(cs11) & cs11 >= 12,
         slow13 = !is.na(cs13) & cs13 >= 12,
         p11 = (low11_19 == 1) | slow11,
         p13 = (low13_19 == 1) | slow13)
cc <- d %>% filter(!is.na(cs11) & !is.na(cs13))
cat(sprintf("complete chair-stand n=%d | labels 2011=%d | not reproduced=%d (%.2f%%) | reproduced=%.2f%%\n",
            nrow(cc), sum(cc$p11), sum(cc$p11 & !cc$p13),
            100*mean(!cc$p13[cc$p11]), 100*mean(cc$p13[cc$p11])))
cat(sprintf("all 6052: labels 2011=%d | not reproduced=%d (%.2f%%) | reproduced=%.2f%%\n",
            sum(d$p11), sum(d$p11 & !d$p13),
            100*mean(!d$p13[d$p11]), 100*mean(d$p13[d$p11])))
write.csv(data.frame(
  sample=c("all 6052","chair-stand complete case"),
  n=c(nrow(d),nrow(cc)),
  labels2011=c(sum(d$p11),sum(cc$p11)),
  not_reproduced=c(sum(d$p11&!d$p13),sum(cc$p11&!cc$p13)),
  reversion_pct=c(100*mean(!d$p13[d$p11]),100*mean(!cc$p13[cc$p11]))),
  file.path("out", "JFA_full_awgs_completecase.csv"), row.names=FALSE)
