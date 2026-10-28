# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# Descriptive comparison behind Section 4.5: the participants excluded for a missing 2013
# ADL block against those retained. Section 1 counts every participant with a grip pattern;
# section 2 restricts to those with age and sex, which is the 2,612 that can be re-added.
## ---- 1. all with a grip pattern -------------------------------------------
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(!is.na(grp19), !is.na(age13), !is.na(male))
b$grp <- relevel(factor(b$grp19, levels=GRP), ref="never")
b$y   <- factor(b$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
b$adl_missing <- is.na(b$a13)
b$adl_disabled <- !is.na(b$a13) & b$a13==1
cat("== all with grip pattern + age + sex:", nrow(b), "\n")
print(b %>% group_by(adl_missing) %>%
  summarise(n=n(), age=mean(age13), male=100*mean(male), grip11=mean(g11,na.rm=TRUE),
            grip13=mean(g13,na.rm=TRUE), low11=100*mean(low11_19,na.rm=TRUE),
            outcome=sum(!is.na(oc)), deaths=sum(y=="dead",na.rm=TRUE),
            disabled=sum(y=="disabled",na.rm=TRUE), .groups="drop"))
cat("\n== the 2,617: what they contribute\n")
x <- b %>% filter(adl_missing)
cat("n =", nrow(x), "| with 2018 outcome =", sum(!is.na(x$oc)),
    "| died =", sum(x$y=="dead",na.rm=TRUE), "| disabled =", sum(x$y=="disabled",na.rm=TRUE),
    "| never-low =", sum(x$grp=="never"), "\n")
cat("\n== retained (a13==0)\n")
r <- b %>% filter(!is.na(a13) & a13==0)
cat("n =", nrow(r), "| age", round(mean(r$age13),1), "| male%", round(100*mean(r$male),1),
    "| grip11", round(mean(r$g11,na.rm=TRUE),1), "| low11%", round(100*mean(r$low11_19,na.rm=TRUE),1), "\n")

## ---- 2. the 2,612 with age and sex ---------------------------------------
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
# the Figure-1 exclusion: grip at both waves, age+sex available, ADL block entirely missing
x <- m %>% filter(!is.na(g11) & !is.na(g13), !is.na(age13), !is.na(male), is.na(a13))
cat("n =", nrow(x), "\n")
cat(sprintf("age %.1f | male %.1f%% | grip11 %.1f kg | low11 %.1f%%\n",
  mean(x$age13), 100*mean(x$male), mean(x$g11,na.rm=TRUE), 100*mean(x$low11_19,na.rm=TRUE)))
r <- m %>% filter(!is.na(g11) & !is.na(g13), !is.na(age13), !is.na(male), a13 %in% 0)
cat(sprintf("retained: age %.1f | male %.1f%% | grip11 %.1f kg | low11 %.1f%% (n=%d)\n",
  mean(r$age13), 100*mean(r$male), mean(r$g11,na.rm=TRUE), 100*mean(r$low11_19,na.rm=TRUE), nrow(r)))