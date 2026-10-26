suppressMessages({library(dplyr)})
m <- readRDS("out/jamda_full.rds")
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
nv <- base %>% filter(grp19 == "never")
set.seed(20261024)
sim1 <- function(){
  res <- list()
  for(sx in c(1,0)){
    nvx <- nv[nv$male==sx,]; sde <- sd(nvx$g13-nvx$g11)/sqrt(2)
    d   <- base[base$male==sx,]; mu <- mean(d$g11); sdo <- sd(d$g11)
    r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if(sx==1) 28 else 18
    Tst <- mu + r*(d$g11-mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r),1e-6)))
    o13 <- Tst + rnorm(nrow(d), 0, sde)
    o13 <- (o13-mean(o13))/sd(o13)*sd(d$g13) + mean(d$g13)   # match mean AND sd
    l1 <- d$g11 < cut; l3 <- o13 < cut
    res[[length(res)+1]] <- c(rev=1-mean(l3[l1]), inc=mean(l3[!l1]))
  }
  colMeans(do.call(rbind,res))
}
R <- replicate(400, sim1())
obs_rev <- 1-mean(base$low13_19[base$low11_19==1]==1)
obs_inc <- mean(base$low13_19[base$low11_19==0]==1)
cat(sprintf("observed:  reversion %.3f | incidence %.3f\n", obs_rev, obs_inc))
cat(sprintf("null:      reversion %.3f (%.3f-%.3f) | incidence %.3f (%.3f-%.3f)\n",
    mean(R["rev",]), quantile(R["rev",],.025), quantile(R["rev",],.975),
    mean(R["inc",]), quantile(R["inc",],.025), quantile(R["inc",],.975)))
cat(sprintf("excess reversion beyond noise: %+.1f pp | observed/null ratio %.2f\n",
    100*(obs_rev-mean(R["rev",])), obs_rev/mean(R["rev",])))
cat(sprintf("incidence ratio observed/null: %.2f\n", obs_inc/mean(R["inc",])))
