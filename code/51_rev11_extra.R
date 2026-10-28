# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
#
# ICC(2,1) and Pearson r for the never-low group in both intervals (Supplementary
# Table S7 Panel B).
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); library(dplyr); 
m <- readRDS("out/jamda_full.rds"); sink("out/51_rev11_extra.txt", split=TRUE)
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
nv <- base %>% filter(grp19=="never")
# ICC(2,1) from a subjects x occasions ANOVA without replication.
# Computed from the sums of squares directly: aov() builds a design matrix with one column per
# subject, which does not fit in memory at this n, and the two routes give identical values.
icc21 <- function(y1, y2){
  ok <- !is.na(y1) & !is.na(y2); y1 <- y1[ok]; y2 <- y2[ok]; n <- length(y1)
  y <- c(y1, y2); gm <- mean(y)
  ss_s <- 2*sum((tapply(y, rep(seq_len(n), 2), mean) - gm)^2)
  ss_o <- n*sum((tapply(y, rep(1:2, each = n), mean) - gm)^2)
  ss_e <- sum((y - gm)^2) - ss_s - ss_o
  msb <- ss_s/(n-1); mse <- ss_e/(n-1)
  c(n = n, msb = msb, mse = mse, icc = (msb-mse)/(msb+mse))
}
r1 <- icc21(nv$g11, nv$g13)
cat(sprintf("never-low n=%d | MSB=%.2f MSE=%.2f | ICC(2,1)=%.3f | ICC(1,1)=%.3f\n",
            r1[["n"]], r1[["msb"]], r1[["mse"]], r1[["icc"]], r1[["icc"]]))
cat(sprintf("Pearson r = %.3f\n", cor(nv$g11, nv$g13)))
# 2013-2015 ICC
b2 <- base[!is.na(base$g15), ]; b2$low15 <- as.integer(b2$g15 < ifelse(b2$male==1,28,18))
nv2 <- b2 %>% filter(low13_19==0, low15==0)
r2 <- icc21(nv2$g13, nv2$g15)
cat(sprintf("2013-2015 never-low n=%d | ICC(2,1)=%.3f | r=%.3f\n", r2[["n"]], r2[["icc"]], cor(nv2$g13,nv2$g15)))
# selection at 2013 on the 9,125 base
two <- m %>% filter(!is.na(g11), !is.na(g13))
cat("\ntwo-wave measured n =", nrow(two), " | a13:", paste(names(table(two$a13, useNA="ifany")), table(two$a13,useNA="ifany"), collapse=" / "), "\n")
two$excl <- as.integer(!(two$a13 %in% 0))
two$low11p <- 100*two$low11_19
for (v in c("age13","g13","g11","bmi","nchronic","cesd","low11p")) {
  cat(sprintf("  %-14s included %.1f | excluded %.1f | p=%.3g\n", v,
      mean(two[[v]][two$excl==0],na.rm=TRUE), mean(two[[v]][two$excl==1],na.rm=TRUE),
      tryCatch(t.test(two[[v]][two$excl==0], two[[v]][two$excl==1])$p.value, error=function(e) NA)))
}
cat(sprintf("  male %%        included %.1f | excluded %.1f | p=%.3g\n",
    100*mean(two$male[two$excl==0]), 100*mean(two$male[two$excl==1]),
    prop.test(c(sum(two$male[two$excl==0]),sum(two$male[two$excl==1])),
              c(sum(two$excl==0),sum(two$excl==1)))$p.value))
cat(sprintf("  died by 2018  included %.1f | excluded %.1f\n",
    100*mean(two$dead18[two$excl==0],na.rm=TRUE), 100*mean(two$dead18[two$excl==1],na.rm=TRUE)))
# 2013-2015 agreement numbers, formatted
agree <- function(la, lb){ ok<-!is.na(la)&!is.na(lb); a<-sum(la[ok]==1&lb[ok]==1); b<-sum(la[ok]==1&lb[ok]==0)
  cc<-sum(la[ok]==0&lb[ok]==1); d<-sum(la[ok]==0&lb[ok]==0); n<-a+b+cc+d
  po<-(a+d)/n; pe<-((a+b)*(a+cc)+(cc+d)*(b+d))/n^2
  c(n=n,a=a,b=b,c=cc,d=d,po=100*po,pe=100*pe,k=(po-pe)/(1-pe),pabak=100*(2*po-1),
    pa=100*2*a/(2*a+b+cc), na=100*2*d/(2*d+b+cc), pi=abs(a-d)/n, bi=abs(b-cc)/n,
    rev=100*b/(a+b), inc=100*cc/(cc+d)) }
print(round(agree(b2$low13_19,b2$low15),3))
print(round(agree(base$low11_19,base$low13_19),3))
# weighted? unweighted only. weighted version for 2013-2015
sink(); cat("ok\n")
