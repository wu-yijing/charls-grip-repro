suppressMessages({library(haven); library(dplyr)})
m <- readRDS("out/jamda_full.rds")
cat("=== 1. 表 2 各组事件数（2019 切点，2018 结局可观测者）===\n")
b <- m %>% filter(a13==0, !is.na(grp19), !is.na(age13), !is.na(male), !is.na(oc))
s <- b %>% group_by(grp19) %>% summarise(n=n(), alive=sum(oc==0), disabled=sum(oc==1), dead=sum(oc==2))
print(as.data.frame(s), row.names=FALSE)
cat("合计: n=", nrow(b), " alive=", sum(b$oc==0), " disabled=", sum(b$oc==1), " dead=", sum(b$oc==2), "\n")
cat("各组事件总数 disabled =", sum(b$oc==1), " dead =", sum(b$oc==2), "\n\n")

cat("=== 2. CES-D 各条目均值（检验是否已反向编码）===\n")
CESD <- paste0("dc01", 1:8)
h <- read_dta("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", col_select = all_of(c("ID", CESD)))
for (v in CESD) {
  lab <- attr(h[[v]], "label"); x <- as.numeric(zap_labels(h[[v]]))
  cat(sprintf("%-6s mean=%.2f  pct4=%.1f%%  %s\n", v, mean(x, na.rm=TRUE),
      100*mean(x==4, na.rm=TRUE), lab))
}
