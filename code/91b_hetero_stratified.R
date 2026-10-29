# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male)) %>%
  mutate(cut=ifelse(male==1,28,18), d=g11-cut, low11=as.numeric(g11<cut), low13=as.numeric(g13<cut))
brks <- c(-Inf,-12,-7,-4,-2,0,2,4,7,12,Inf); b$str <- cut(b$d, brks, right=FALSE); lev<-levels(b$str)
sd_chg <- sapply(lev, function(s) sd(b$g13[b$str==s]-b$g11[b$str==s]))
n_s <- sapply(lev, function(s) sum(b$str==s))
mu <- mean(b$g11); varX <- var(b$g11)
r_paper <- 0.674; se_paper <- sqrt(1-r_paper)*sqrt(varX)   # 5.54
base_se <- sd_chg/sqrt(2)
k <- se_paper / (sum(n_s*base_se)/sum(n_s))
se_s <- base_se*k
sT2 <- varX - sum(n_s*se_s^2)/sum(n_s)
r_s <- sT2/(sT2+se_s^2)
cat(sprintf("se_paper=%.3f k=%.3f sigma_T=%.2f\n", se_paper, k, sqrt(sT2)))
print(data.frame(stratum=lev, n=as.integer(n_s), sd_change=round(as.numeric(sd_chg),2),
                 sigma_e=round(as.numeric(se_s),2), r=round(as.numeric(r_s),3)))
sim <- function(se,r,target,reps=400,seed=20261028){
  set.seed(seed); n<-nrow(b); out<-matrix(0,reps,2)
  for(kk in 1:reps){
    lat<-numeric(n); g<-numeric(n)
    for(s in seq_along(lev)){
      idx<-which(b$str==lev[s])
      lat[idx]<-mu+r[s]*(b$g11[idx]-mu)+rnorm(length(idx),0,sqrt(r[s]*(1-r[s]))*sqrt(varX))
      g[idx]<-lat[idx]+rnorm(length(idx),0,se[s])
    }
    f<-function(dl) mean(as.numeric((g+dl)<b$cut))-target
    dl<-uniroot(f,c(-30,30))$root; gg<-g+dl
    out[kk,1]<-mean(gg[b$low11==1]>=b$cut[b$low11==1]); out[kk,2]<-mean(gg[b$low11==0]<b$cut[b$low11==0])
  }
  c(colMeans(out), quantile(out[,1],c(.025,.975)))
}
obs<-mean(b$g13[b$low11==1]>=b$cut[b$low11==1])
# paper spec: homogeneous r=0.674, se=5.54
for(tg in c("2011","2013")){
  t0<-if(tg=="2011") mean(b$low11) else mean(b$low13)
  s<-sim(rep(se_paper,10), rep(r_paper,10), t0)
  cat(sprintf("HOMO(paper) %s: rev %.4f (%.4f-%.4f) %.0f%% of obs; inc %.4f\n",tg,s[1],s[3],s[4],100*s[1]/obs,s[2]))
}
for(tg in c("2011","2013")){
  t0<-if(tg=="2011") mean(b$low11) else mean(b$low13)
  s<-sim(se_s, r_s, t0)
  cat(sprintf("HETERO      %s: rev %.4f (%.4f-%.4f) %.0f%% of obs; inc %.4f\n",tg,s[1],s[3],s[4],100*s[1]/obs,s[2]))
}
