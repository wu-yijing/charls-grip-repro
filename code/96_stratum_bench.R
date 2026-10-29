# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
# 96_stratum_bench.R -- 主规格下的分层检验
# 规格与 31_rev9_errorsim.R V1 完全一致：误差 SD 取自 never-low（按性别），
# 矩（mu, sdo）取自全样本（按性别），只把被模拟的人限制在层内。
# 校准：g11 = 2011 边际（正文主校准，全样本应复现 59.5%）；g13 = 2013 边际（55.6%）
suppressMessages(library(dplyr))
set.seed(20261024)

m <- readRDS("out/jamda_full.rds")
base <- m %>% filter(a13 == 0, !is.na(grp19), !is.na(age13), !is.na(male))

sim <- function(dat, calib = c("g11","g13"), R = 400) {
  calib <- match.arg(calib)
  rev <- numeric(0)
  for (sx in c(1, 0)) {
    d  <- dat[dat$male == sx, ]
    nv <- base[base$grp19 == "never" & base$male == sx, ]          # 全样本 never-low
    if (nrow(d) < 30) next
    sde <- sd(nv$g13 - nv$g11)/sqrt(2)
    mu  <- mean(base$g11[base$male == sx]); sdo <- sd(base$g11[base$male == sx])
    r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if (sx == 1) 28 else 18
    ref <- if (calib == "g11") d$g11 else d$g13
    l1  <- d$g11 < cut
    if (!any(l1)) next
    for (b in seq_len(R)) {
      Tst <- mu + r*(d$g11 - mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r), 1e-6)))
      o13 <- Tst + rnorm(nrow(d), 0, sde)
      o13 <- (o13 - mean(o13))/sd(o13)*sd(ref) + mean(ref)
      rev <- c(rev, 1 - mean((o13 < cut)[l1]))
    }
  }
  rev
}

for (lab in c("all","<65",">=65")) {
  dat <- switch(lab, "all"=base, "<65"=base[base$age13<65,], ">=65"=base[base$age13>=65,])
  l1 <- dat$low11_19 == 1
  obs <- 1 - mean(dat$low13_19[l1])
  e11 <- sim(dat,"g11"); e13 <- sim(dat,"g13")
  cat(sprintf("%-4s n=%4d low11=%3d | observed %5.1f%% | g11 exp %5.1f%% (%4.1f-%4.1f) ratio %.2f | g13 exp %5.1f%% ratio %.2f\n",
    lab, nrow(dat), sum(l1), 100*obs,
    100*mean(e11), 100*quantile(e11,.025), 100*quantile(e11,.975), obs/mean(e11),
    100*mean(e13), obs/mean(e13)))
}
