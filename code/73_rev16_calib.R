# 73_rev16_calib.R -- calibration of the no-change model: matching to the 2013 or the 2011 marginal
#
# Supplementary Table S14. The three rows differ only in the final moment match applied to the
# simulated 2013 value. The error term is the never-low SEM (per sex) and the reliability is
# 1 - SEM^2 / SD(g11)^2 within sex, i.e. 0.66 in men and 0.51 in women (Table S32, Panel A);
# there is no single pooled 0.674 in this script.
#
#   C1  mean and SD of the OBSERVED 2013 value   (primary; imposes the observed 0.55 kg decline)
#   C2  mean and SD of the 2011 value            (strict no change; reproduces 59.5% by construction)
#   C3  mean of 2011, SD of observed 2013        (mean-only correction)
#
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages(library(dplyr))
set.seed(20261024)

m <- readRDS("out/jamda_full.rds")
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
nv   <- base %>% filter(grp19 == "never")
cat("n =", nrow(base), "| never-low =", nrow(nv), "\n")
cat(sprintf("mean 2011 = %.2f kg | mean 2013 = %.2f kg | difference = %.2f kg\n",
            mean(base$g11), mean(base$g13), mean(base$g13) - mean(base$g11)))

sim1 <- function(cal){
  res <- list()
  for (sx in c(1,0)) {
    nvx <- nv[nv$male==sx, ]
    sde <- sd(nvx$g13 - nvx$g11)/sqrt(2)
    d   <- base[base$male==sx, ]
    mu  <- mean(d$g11); sdo <- sd(d$g11)
    r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if (sx==1) 28 else 18
    Tst <- mu + r*(d$g11-mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r), 1e-6)))
    o13 <- Tst + rnorm(nrow(d), 0, sde)
    o13 <- switch(cal,
      C1 = (o13 - mean(o13))/sd(o13)*sd(d$g13) + mean(d$g13),
      C2 = (o13 - mean(o13))/sd(o13)*sd(d$g11) + mean(d$g11),
      C3 = (o13 - mean(o13))/sd(o13)*sd(d$g13) + mean(d$g11))
    l1 <- d$g11 < cut; l3 <- o13 < cut
    res[[length(res)+1]] <- c(rev = 1-mean(l3[l1]), inc = mean(l3[!l1]))
  }
  colMeans(do.call(rbind, res))
}

obs_rev <- 1 - mean(base$low13_19[base$low11_19==1] == 1)
obs_inc <- mean(base$low13_19[base$low11_19==0] == 1)
cat(sprintf("\nobserved: reversion %.4f | incidence %.4f\n\n", obs_rev, obs_inc))

out <- list()
for (cal in c("C1","C2","C3")) {
  R <- replicate(400, sim1(cal))
  er <- mean(R["rev",]); ei <- mean(R["inc",])
  out[[cal]] <- data.frame(
    calibration = cal,
    exp_rev = er, rev_lo = quantile(R["rev",], .025), rev_hi = quantile(R["rev",], .975),
    rev_ratio = obs_rev/er, rev_excess_pp = 100*(obs_rev-er),
    exp_inc = ei, inc_lo = quantile(R["inc",], .025), inc_hi = quantile(R["inc",], .975),
    inc_ratio = obs_inc/ei)
  cat(sprintf("%s  expected reversion %.1f%% (%.1f-%.1f) | obs/exp %.2f | excess %+.1f pp\n",
      cal, 100*er, 100*quantile(R["rev",],.025), 100*quantile(R["rev",],.975), obs_rev/er,
      100*(obs_rev-er)))
  cat(sprintf("    expected incidence %.1f%% (%.1f-%.1f) | obs/exp %.2f\n",
      100*ei, 100*quantile(R["inc",],.025), 100*quantile(R["inc",],.975), obs_inc/ei))
}
res <- do.call(rbind, out)
write.csv(res, "out/73_calib.csv", row.names=FALSE)
sink("out/73_calib.txt")
cat("S-2 calibration sensitivity of the no-change model\n")
cat(sprintf("n = %d; observed reversion %.3f; observed incidence %.3f; mean 2011 %.2f; mean 2013 %.2f\n\n",
            nrow(base), obs_rev, obs_inc, mean(base$g11), mean(base$g13)))
print(res, digits = 3, row.names = FALSE)
sink()
cat("\nwritten out/73_calib.csv / .txt\n")
