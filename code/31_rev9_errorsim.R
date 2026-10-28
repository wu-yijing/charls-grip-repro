# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
#
# Three specifications of the error term in the no-change benchmark (Supplementary
# Table S9): constant SD from the never-low group (published), constant SD anchored to
# the persistent group, and an SD proportional to the 2011 value.
# 31_rev9_errorsim.R -- rev9 addition for reviewer M3-1
# The published null simulation takes the error SD from the NEVER-LOW group and
# applies it unchanged to everyone. Here the same generative model is re-run with:
#   V1  primary  : constant SD estimated in never-low  (as published)
#   V2  low-anchor: constant SD estimated in persistent (both waves low)
#   V3  proportional: SD_i = CV * own 2011 value, CV from never-low
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
library(dplyr)

m <- readRDS("out/jamda_full.rds")
set.seed(20261024)

base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
nv   <- base %>% filter(grp19 == "never")
ps   <- base %>% filter(grp19 == "persistent")
cat("n base =", nrow(base), "| never =", nrow(nv), "| persistent =", nrow(ps), "\n")

sem_nv <- mean(c(sd(nv$g13[nv$male==1]-nv$g11[nv$male==1]),
                 sd(nv$g13[nv$male==0]-nv$g11[nv$male==0])))/sqrt(2)
sem_ps <- mean(c(sd(ps$g13[ps$male==1]-ps$g11[ps$male==1]),
                 sd(ps$g13[ps$male==0]-ps$g11[ps$male==0])))/sqrt(2)
cv     <- sem_nv/mean(nv$g11)
cat(sprintf("SEM never-low = %.2f kg | SEM persistent = %.2f kg | CV = %.3f\n", sem_nv, sem_ps, cv))

sim1 <- function(mode){
  res <- list()
  for (sx in c(1,0)) {
    nvx <- nv[nv$male==sx, ]
    sde <- sd(nvx$g13-nvx$g11)/sqrt(2)
    d   <- base[base$male==sx, ]
    mu  <- mean(d$g11); sdo <- sd(d$g11)
    r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if (sx==1) 28 else 18
    Tst <- mu + r*(d$g11-mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r),1e-6)))
    e   <- switch(mode,
      V1 = rnorm(nrow(d), 0, sde),
      V2 = rnorm(nrow(d), 0, sem_ps),
      V3 = rnorm(nrow(d), 0, cv*d$g11))
    o13 <- Tst + e
    o13 <- (o13-mean(o13))/sd(o13)*sd(d$g13) + mean(d$g13)
    l1 <- d$g11 < cut; l3 <- o13 < cut
    res[[length(res)+1]] <- c(rev=1-mean(l3[l1]), inc=mean(l3[!l1]))
  }
  colMeans(do.call(rbind, res))
}

obs_rev <- 1-mean(base$low13_19[base$low11_19==1]==1)
obs_inc <- mean(base$low13_19[base$low11_19==0]==1)
cat(sprintf("observed: reversion %.4f | incidence %.4f\n\n", obs_rev, obs_inc))

for (mode in c("V1","V2","V3")) {
  R <- replicate(400, sim1(mode))
  er <- mean(R["rev",]); ei <- mean(R["inc",])
  cat(sprintf("%s  expected reversion %.1f%% (%.1f-%.1f) | ratio obs/exp %.2f | implied noise share %.0f%%\n",
      mode, 100*er, 100*quantile(R["rev",],.025), 100*quantile(R["rev",],.975),
      obs_rev/er, 100*min(1, obs_rev/er)))
  cat(sprintf("     expected incidence %.1f%% (%.1f-%.1f) | ratio obs/exp %.2f\n",
      100*ei, 100*quantile(R["inc",],.025), 100*quantile(R["inc",],.975), obs_inc/ei))
}
