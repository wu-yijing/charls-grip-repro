suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
for (v in c("19","25")) {
  b <- m %>% filter(a13==0, !is.na(.data[[paste0("grp",v)]]), !is.na(age13), !is.na(male))
  t <- table(factor(b[[paste0("low11_",v)]],0:1), factor(b[[paste0("low13_",v)]],0:1))
  po <- sum(diag(t))/sum(t); pe <- sum(rowSums(t)*colSums(t))/sum(t)^2
  cat(sprintf("cutoffs %s: n=%d cell=[%d %d; %d %d] Po=%.3f Pe=%.3f kappa=%.3f\n",
      v, sum(t), t[1,1], t[1,2], t[2,1], t[2,2], po, pe, (po-pe)/(1-pe)))
}
