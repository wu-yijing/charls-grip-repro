suppressMessages({library(dplyr)})
m <- readRDS("out/jamda_full.rds")
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
nv <- base %>% filter(grp19 == "never")
obs_rev <- 1 - mean(base$low13_19[base$low11_19 == 1] == 1)

sim <- function(fac, reps = 300){
  set.seed(20261024)
  R <- replicate(reps, {
    res <- list()
    for (sx in c(1, 0)) {
      nvx <- nv[nv$male == sx, ]
      sde <- fac * sd(nvx$g13 - nvx$g11) / sqrt(2)
      d   <- base[base$male == sx, ]; mu <- mean(d$g11); sdo <- sd(d$g11)
      r   <- max(0.05, 1 - sde^2 / sdo^2); cut <- if (sx == 1) 28 else 18
      Tst <- mu + r * (d$g11 - mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2 * (1 - r), 1e-6)))
      o13 <- Tst + rnorm(nrow(d), 0, sde)
      o13 <- (o13 - mean(o13)) / sd(o13) * sd(d$g13) + mean(d$g13)
      l1 <- d$g11 < cut; l3 <- o13 < cut
      res[[length(res) + 1]] <- 1 - mean(l3[l1])
    }
    mean(unlist(res))
  })
  c(lo = quantile(R, .025), mean = mean(R), hi = quantile(R, .975))
}
cat("observed reversion:", round(obs_rev, 4), "\n\n")
for (f in c(0.7, 0.8, 0.9, 1.0, 1.2, 1.5)) {
  s <- sim(f)
  cat(sprintf("SEM multiplier %.1f -> null reversion %.3f (%.3f-%.3f) | noise share %.1f%%\n",
              f, s["mean"], s["lo"], s["hi"], 100 * s["mean"] / obs_rev))
}
