# 21_enhanced_A.R -- reproducibility / measurement-error analyses for the Option-A rewrite (v2)
suppressMessages({library(nnet); library(dplyr); library(pROC); library(sandwich); library(MASS)})
m <- readRDS("out/jamda_full.rds")
sink("out/21_enhanced_A_v2.txt", split=TRUE)

cat("#################### A. FLOW ####################\n")
for(nm in c("low-based","grip-based")){
  if(nm=="low-based") cl <- !is.na(m$low11_19) & !is.na(m$low13_19)
  else                cl <- !is.na(m$g11) & !is.na(m$g13)
  cat(nm, "classifiable:", sum(cl),
      "| +age&sex:", sum(cl & !is.na(m$age13) & !is.na(m$male)), "\n")
}
cl <- !is.na(m$low11_19) & !is.na(m$low13_19)
cat("ADL 2013 within low-based classifiable:\n"); print(table(m$a13[cl], useNA="ifany"))
cat("ADL 2013 within classifiable + age/sex:\n")
print(table(m$a13[cl & !is.na(m$age13) & !is.na(m$male)], useNA="ifany"))

base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
cat("analytic:", nrow(base), "| groups:\n"); print(table(base$grp19)[c("never","remitted","incident","persistent")])

cat("\n#################### B. AGREEMENT ####################\n")
agree <- function(l1,l3,tag,w=NULL){
  tb <- table(factor(l1,0:1), factor(l3,0:1)); n <- sum(tb)
  po <- sum(diag(tb))/n; pe <- sum(rowSums(tb)*colSums(tb))/n^2; k <- (po-pe)/(1-pe)
  set.seed(20261024); ks <- numeric(1000); cell <- l1*2 + l3 + 1L
  for(b in 1:1000){ s <- sample.int(n,n,TRUE); t2 <- matrix(tabulate(cell[s],nbins=4),2,2)
    p2 <- sum(diag(t2))/n; e2 <- sum(rowSums(t2)*colSums(t2))/n^2; ks[b] <- (p2-e2)/(1-e2) }
  a<-tb[2,2]; b<-tb[2,1]; c<-tb[1,2]; d<-tb[1,1]
  cat("\n--",tag,"(n =",n,")--\n"); print(tb)
  cat(sprintf("obs %.4f | exp %.4f | kappa %.4f (%.4f-%.4f)\n",po,pe,k,quantile(ks,.025),quantile(ks,.975)))
  cat(sprintf("PABAK %.4f | pos %.4f | neg %.4f | prev idx %.4f | bias idx %.4f\n",
      2*po-1, 2*a/(2*a+b+c), 2*d/(2*d+b+c), abs((a+b)-(a+c))/n, abs(b-c)/n))
}
i_full <- which(cl & !is.na(m$age13) & !is.na(m$male))
agree(as.integer(m$low11_19)[i_full], as.integer(m$low13_19)[i_full], "two-wave classifiable")
i_an <- which(m$a13 %in% 0 & !is.na(m$grp19) & !is.na(m$age13) & !is.na(m$male))
agree(as.integer(m$low11_19)[i_an], as.integer(m$low13_19)[i_an], "analytic")

cat("\n#################### C. CONTINUOUS NOISE (analytic sample) ####################\n")
nv <- base %>% filter(grp19=="never")
for(k in c("unweighted","weighted")){
  w <- if(k=="weighted") ifelse(is.na(nv$w13),0,nv$w13) else rep(1,nrow(nv))
  ok<- w>0
  mu1 <- weighted.mean(nv$g11[ok], w[ok]); mu3 <- weighted.mean(nv$g13[ok], w[ok])
  s1 <- sqrt(sum(w[ok]*(nv$g11[ok]-mu1)^2)/sum(w[ok])); s3 <- sqrt(sum(w[ok]*(nv$g13[ok]-mu3)^2)/sum(w[ok]))
  dd <- nv$g13-nv$g11; md <- weighted.mean(dd[ok],w[ok]); sd_ <- sqrt(sum(w[ok]*(dd[ok]-md)^2)/sum(w[ok]))
  r <- sum(w[ok]*(nv$g11[ok]-mu1)*(nv$g13[ok]-mu3))/sqrt(sum(w[ok]*(nv$g11[ok]-mu1)^2)*sum(w[ok]*(nv$g13[ok]-mu3)^2))
  cat(sprintf("\n%s never-low (n=%d): 2011 %.2f (%.2f) | 2013 %.2f (%.2f) | change %+.2f (%.2f) | r %.3f\n",
      k, sum(ok), mu1,s1,mu3,s3,md,sd_,r))
  cat(sprintf("  SEM %.2f | MDC95 %.2f | LoA %.2f to %.2f | outside %.1f%%\n",
      sd_/sqrt(2), 1.96*sd_, md-1.96*sd_, md+1.96*sd_, 100*mean(abs(dd-md)>1.96*sd_)))
}
dd <- base$g13-base$g11
cat(sprintf("\nanalytic sample change: mean %+.2f SD %.2f\n", mean(dd), sd(dd)))
cat(sprintf("Bland-Altman slope (diff on mean) %.3f (SE %.3f)\n",
    coef(lm(I(g13-g11)~I((g11+g13)/2),data=base))[2],
    summary(lm(I(g13-g11)~I((g11+g13)/2),data=base))$coef[2,2]))
# by sex
for(s in c(1,0)){ d<-nv[nv$male==s,]; dd2<-d$g13-d$g11
  cat(sprintf("  sex=%s: change %+.2f SD %.2f MDC95 %.2f\n", s, mean(dd2), sd(dd2), 1.96*sd(dd2))) }

cat("\n#################### D. MEASUREMENT-ERROR NULL (analytic sample) ####################\n")
set.seed(20261024)
sim1 <- function(){
  res <- list()
  for(sx in c(1,0)){
    nvx <- nv[nv$male==sx,]; sde <- sd(nvx$g13-nvx$g11)/sqrt(2)
    d   <- base[base$male==sx,]; mu <- mean(d$g11); sdo <- sd(d$g11)
    r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if(sx==1) 28 else 18
    Tst <- mu + r*(d$g11-mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r),1e-6)))
    o13 <- Tst + rnorm(nrow(d), 0, sde)
    o13 <- o13 - mean(o13) + mean(d$g13)      # calibrate to observed wave mean
    l1 <- d$g11 < cut; l3 <- o13 < cut
    res[[length(res)+1]] <- c(rev=mean(l3[l1]), inc=mean(l3[!l1]))
  }
  colMeans(do.call(rbind,res))
}
R <- replicate(400, sim1())
obs_rev <- mean(base$low13_19[base$low11_19==1]==1)
obs_inc <- mean(base$low13_19[base$low11_19==0]==1)
cat(sprintf("observed (analytic): reversion %.3f | incidence %.3f\n", obs_rev, obs_inc))
cat(sprintf("expected under noise-only null: reversion %.3f (2.5-97.5: %.3f-%.3f) | incidence %.3f (%.3f-%.3f)\n",
    mean(R["rev",]), quantile(R["rev",],.025), quantile(R["rev",],.975),
    mean(R["inc",]), quantile(R["inc",],.025), quantile(R["inc",],.975)))
cat(sprintf("observed/expected: reversion %.2f | incidence %.2f\n", obs_rev/mean(R["rev",]), obs_inc/mean(R["inc",])))
cat(sprintf("excess reversion not explained by noise: %.1f pp\n", 100*(obs_rev-mean(R["rev",]))))

cat("\n#################### E. ORDINAL-CONSTRAINT LR TEST ####################\n")
d0 <- base; d0$y <- factor(d0$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
d0$grp <- relevel(factor(d0$grp19, levels=c("never","remitted","incident","persistent")), ref="never")
d0$nr <- d0$low11_19 + d0$low13_19; d0$wsc <- d0$w13/mean(d0$w13, na.rm=TRUE)
COV <- c("age10","male","bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4")
dc <- d0[complete.cases(d0[,c("y",COV,"grp","nr")]),]
m4 <- multinom(as.formula(paste("y ~",paste(c(COV,"grp"),collapse="+"))), data=dc, weights=dc$wsc, trace=FALSE, maxit=500)
m2 <- multinom(as.formula(paste("y ~",paste(c(COV,"nr"), collapse="+"))), data=dc, weights=dc$wsc, trace=FALSE, maxit=500)
ll4 <- as.numeric(logLik(m4)); ll2 <- as.numeric(logLik(m2)); df <- 4
cat(sprintf("4-group %.2f | ordinal %.2f | LR chi2 %.3f df %d P %.4f\n", ll4, ll2, 2*(ll4-ll2), df,
    pchisq(2*(ll4-ll2), df, lower.tail=FALSE)))

cat("\n#################### F. STANDARDISED RISK & RISK DIFFERENCE at 65 ####################\n")
Xpred <- function(d, co){
  mu <- sapply(d[,COV], mean, na.rm=TRUE)
  out <- list()
  for(g in c("never","remitted","incident","persistent")){
    x <- c(1, mu, as.numeric(g==c("remitted","incident","persistent")))
    eta <- matrix(x, nrow=1) %*% t(matrix(co, nrow=2, byrow=TRUE))
    p <- exp(c(0,eta)); p <- p/sum(p); out[[g]] <- p
  }
  do.call(rbind, out)
}
co <- coef(m4); rownames(co) <- c("disabled","dead")
P <- Xpred(dc, as.numeric(t(co))); colnames(P) <- c("alive","disabled","dead")
rownames(P) <- c("never","remitted","incident","persistent")
cat("complete-case weighted standardized risks at age 65 (reference profile):\n"); print(round(P*100,1))
cat(sprintf("disability risk difference persistent-never: %+.1f pp\n", 100*(P["persistent","disabled"]-P["never","disabled"])))
cat(sprintf("mortality  risk difference persistent-never: %+.1f pp\n", 100*(P["persistent","dead"]-P["never","dead"])))
cat(sprintf("disability risk difference remitted-never:  %+.1f pp\n", 100*(P["remitted","disabled"]-P["never","disabled"])))
cat(sprintf("mortality  risk difference remitted-never:  %+.1f pp\n", 100*(P["remitted","dead"]-P["never","dead"])))

cat("\n#################### G. DECISION CURVE (disability) ####################\n")
dd2 <- dc[dc$y!="dead",]; dd2$dis <- as.numeric(dd2$y=="disabled")
fA <- as.formula(paste("dis ~", paste(c(COV,"g13"), collapse="+")))
fB <- as.formula(paste("dis ~", paste(c(COV,"g13","grp"), collapse="+")))
gA <- glm(fA, data=dd2, family=binomial); gB <- glm(fB, data=dd2, family=binomial)
pA <- fitted(gA); pB <- fitted(gB); y <- dd2$dis; n <- length(y)
nb <- function(p,pt){ tp<-mean(p>=pt & y==1); fp<-mean(p>=pt & y==0); tp - fp*pt/(1-pt) }
ths <- c(0.05,0.10,0.15,0.20,0.25,0.30)
cat(sprintf("n=%d events=%d\n", n, sum(y)))
cat(sprintf("%-8s %10s %10s %10s %10s\n","thresh","treat-all","modelA","modelB2","B2-A"))
for(pt in ths){ cat(sprintf("%-8.2f %10.4f %10.4f %10.4f %+10.4f\n", pt,
  nb(rep(1,n),pt), nb(pA,pt), nb(pB,pt), nb(pB,pt)-nb(pA,pt))) }
cat(sprintf("max net benefit gain from adding the pattern: %+.4f\n",
    max(sapply(seq(0.03,0.40,0.01), function(pt) nb(pB,pt)-nb(pA,pt)))))

cat("\n#################### H. RETEST YIELD / 2015 ####################\n")
cat(sprintf("2011 low readings (two-wave classifiable): %d | persistent: %d | retests per case: %.2f\n",
    sum(m$low11_19[i_full]==1, na.rm=TRUE), sum(m$low11_19[i_full]==1 & m$low13_19[i_full]==1, na.rm=TRUE),
    sum(m$low11_19[i_full]==1,na.rm=TRUE)/sum(m$low11_19[i_full]==1 & m$low13_19[i_full]==1,na.rm=TRUE)))
p <- base %>% filter(grp19=="persistent", !is.na(g15))
lo15 <- p$g15 < ifelse(p$male==1,28,18)
cat(sprintf("persistent with 2015 reading: %d | still low: %d (%.1f%%) | normal: %d (%.1f%%)\n",
    nrow(p), sum(lo15), 100*mean(lo15), sum(!lo15), 100*mean(!lo15)))
# among all 2011-low with 2015 data: how many low at all three waves
q <- base %>% filter(low11_19==1, !is.na(g15)); lo15b <- q$g15 < ifelse(q$male==1,28,18)
cat(sprintf("2011-low with 2015 reading: %d | low at 3/3: %d (%.1f%%) | low at 1/3: %d (%.1f%%)\n",
    nrow(q), sum(q$low13_19==1 & lo15b), 100*mean(q$low13_19==1 & lo15b),
    sum(q$low13_19==0 & !lo15b), 100*mean(q$low13_19==0 & !lo15b)))
sink(); cat("done\n")
