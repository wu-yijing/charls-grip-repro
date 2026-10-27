# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 50_rev11_review.R -- analyses added in response to the JNHA peer review (rev11)
#  (A) 2013-2015 replication of the two-year agreement (review M1b)
#  (B) mean-of-two-readings model for incremental value (review M2)
#  (C) coefficient of the current grip value in model A (review M8)
#  (D) included vs excluded comparison for missing 2013 ADL (review M5)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
set.seed(20261024)
sink("out/50_rev11_review.txt", split=TRUE)

cat("=============================================================\n")
cat("(0) availability of the 2015 grip\n")
cat("=============================================================\n")
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
cat("analytic sample n =", nrow(base), "\n")
cat("g15 non-missing in analytic sample:", sum(!is.na(base$g15)), "\n")
cat("g15 non-missing and alive at 2015 (oc != 2 or dead18==0):",
    sum(!is.na(base$g15) & (is.na(base$dead18) | base$dead18==0)), "\n")
print(table(g15_avail = !is.na(base$g15), grp = base$grp19))

# ---------------------------------------------------------------
cat("\n=============================================================\n")
cat("(A) two-year agreement: 2011-2013 vs 2013-2015\n")
cat("=============================================================\n")

agree <- function(low_a, low_b, dat) {
  ok <- !is.na(low_a) & !is.na(low_b)
  a <- sum(low_a[ok]==1 & low_b[ok]==1); b <- sum(low_a[ok]==1 & low_b[ok]==0)
  cc<- sum(low_a[ok]==0 & low_b[ok]==1); d <- sum(low_a[ok]==0 & low_b[ok]==0)
  n <- a+b+cc+d
  po <- (a+d)/n
  pe <- ((a+b)*(a+cc) + (cc+d)*(b+d))/n^2
  kap<- (po-pe)/(1-pe)
  pa <- 2*a/(2*a+b+cc); na <- 2*d/(2*d+b+cc)
  c(n=n, a=a, b=b, c=cc, d=d, po=po, pe=pe, kappa=kap,
    pabak=2*po-1, pos_agree=pa, neg_agree=na,
    prev_index=abs(a-d)/n, bias_index=abs(b-cc)/n,
    rev=ifelse((a+b)>0, b/(a+b), NA), inc=ifelse((cc+d)>0, cc/(cc+d), NA))
}

r1 <- agree(base$low11_19, base$low13_19, base)
b2 <- base[!is.na(base$g15), ]
b2$low15_19 <- as.integer(b2$g15 < ifelse(b2$male==1, 28, 18))
r2 <- agree(b2$low13_19, b2$low15_19, b2)

fmt <- function(r, lab) {
  cat(sprintf("\n--- %s ---\n", lab))
  cat(sprintf("n = %d | %d/%d; %d/%d; %d/%d; %d/%d\n", r["n"], r["a"], r["a"]+r["b"],
              r["b"], r["a"]+r["b"], r["c"], r["c"]+r["d"], r["d"], r["c"]+r["d"]))
  cat(sprintf("observed agreement %.1f%% | expected %.1f%% | kappa %.3f | PABAK %.3f\n",
              100*r["po"], 100*r["pe"], r["kappa"], r["pabak"]))
  cat(sprintf("positive agreement %.1f%% | negative agreement %.1f%% | prev index %.3f | bias index %.3f\n",
              100*r["pos_agree"], 100*r["neg_agree"], r["prev_index"], r["bias_index"]))
  cat(sprintf("reversion %.1f%% | incidence %.1f%%\n", 100*r["rev"], 100*r["inc"]))
}
fmt(r1, "2011 -> 2013 (analytic sample)")
fmt(r2, "2013 -> 2015 (analytic sample with a 2015 reading)")

# who is in the 2013-2015 analysis vs not
cat("\ncomposition of the 2013-2015 analysis set by 2011-2013 pattern:\n")
print(round(100*prop.table(table(base$grp19)),1))
print(round(100*prop.table(table(b2$grp19)),1))

cat("\nby sex, 2013->2015:\n")
for (sx in c(1,0)) {
  ss <- b2[b2$male==sx, ]
  rr <- agree(ss$low13_19, ss$low15_19, ss)
  cat(sprintf("  %s: n=%d kappa=%.3f po=%.3f rev=%.3f inc=%.3f\n",
              ifelse(sx==1,"men","women"), rr["n"], rr["kappa"], rr["po"], rr["rev"], rr["inc"]))
}

# ---------------------------------------------------------------
cat("\n=============================================================\n")
cat("(A2) no-change benchmark for 2013-2015\n")
cat("=============================================================\n")
nv2 <- b2 %>% filter(low13_19==0, low15_19==0)
cat("never-low 2013-2015 n =", nrow(nv2), "\n")
sd_chg <- sd(nv2$g15 - nv2$g13); sem2 <- sd_chg/sqrt(2)
cat(sprintf("SD(change) = %.2f kg | SEM = %.2f kg | MDC95 = %.2f kg\n",
            sd_chg, sem2, 1.96*sd_chg))
cat(sprintf("between-person SD of 2013 grip in nv2 = %.2f kg | implied reliability %.2f\n",
            sd(nv2$g13), 1 - sem2^2/sd(nv2$g13)^2))
cat(sprintf("mean grip 2013 %.1f | 2015 %.1f | mean change %.1f (SD %.1f)\n",
            mean(nv2$g13), mean(nv2$g15), mean(nv2$g15-nv2$g13), sd(nv2$g15-nv2$g13)))
cat(sprintf("correlation 2013/2015 = %.3f | 2011/2013 = %.3f\n",
            cor(nv2$g13, nv2$g15), cor(m$g11[m$grp19=="never" & m$a13 %in% 0], m$g13[m$grp19=="never" & m$a13 %in% 0])))

sim2 <- function(reps=400) {
  out <- matrix(NA, reps, 2, dimnames=list(NULL, c("rev","inc")))
  for (i in 1:reps) {
    rr <- c()
    for (sx in c(1,0)) {
      nvx <- nv2[nv2$male==sx, ]
      sde <- sd(nvx$g15 - nvx$g13)/sqrt(2)
      d   <- b2[b2$male==sx, ]
      mu  <- mean(d$g13); sdo <- sd(d$g13)
      r   <- max(0.3, 1 - sde^2/sdo^2); cut <- if (sx==1) 28 else 18
      Tst <- mu + r*(d$g13-mu) + rnorm(nrow(d), 0, sqrt(max(sdo^2*(1-r),1e-6)))
      o15 <- Tst + rnorm(nrow(d), 0, sde)
      o15 <- (o15-mean(o15))/sd(o15)*sd(d$g15) + mean(d$g15)
      l3 <- d$low13_19==1; l5 <- o15 < cut
      rr <- rbind(rr, c(1-mean(l5[l3]), mean(l5[!l3])))
    }
    out[i,] <- colMeans(rr)
  }
  out
}
S <- sim2(400)
cat(sprintf("expected reversion %.1f%% (%.1f-%.1f) vs observed %.1f%%; obs/exp %.2f\n",
            100*mean(S[,"rev"]), 100*quantile(S[,"rev"],.025), 100*quantile(S[,"rev"],.975),
            100*r2["rev"], r2["rev"]/mean(S[,"rev"])))
cat(sprintf("expected incidence %.1f%% (%.1f-%.1f) vs observed %.1f%%; obs/exp %.2f\n",
            100*mean(S[,"inc"]), 100*quantile(S[,"inc"],.025), 100*quantile(S[,"inc"],.975),
            100*r2["inc"], r2["inc"]/mean(S[,"inc"])))

# three-wave: what fraction of 2011 lows are low at 2015 too
cat("\nthree-wave persistence (2011 low -> 2013 low -> 2015 low), n with all three:\n")
b3 <- base[!is.na(base$g15), ]
b3$low15_19 <- as.integer(b3$g15 < ifelse(b3$male==1, 28, 18))
print(table(low11=b3$low11_19, low13=b3$low13_19, low15=b3$low15_19))
cat(sprintf("of 2011 lows with all three readings: still low at 2015 %.1f%%\n",
            100*mean(b3$low15_19[b3$low11_19==1])))

# ---------------------------------------------------------------
cat("\n=============================================================\n")
cat("(B) mean-of-two-readings model (review M2) and (C) current-grip OR\n")
cat("=============================================================\n")
b <- get_base(m, "grp19")
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5; b$chg_5 <- (b$g13-b$g11)/5
b$gmean_5 <- (b$g11+b$g13)/10
b$low11 <- b$low11_19; b$low13 <- b$low13_19; b$grp4 <- b$grp
d0 <- b[complete.cases(b[, c(COV,"g13_5","g11_5","gmean_5","chg_5","low11","low13","grp4","oc")]), ]
cat("complete-case analysis set n =", nrow(d0), "\n")

nri_cont <- function(pA,pB,event){up<-pB>pA;dn<-pB<pA;e<-event==1
  (mean(up[e])-mean(dn[e]))+(mean(dn[!e])-mean(up[!e]))}
idi <- function(pA,pB,event){e<-event==1
  (mean(pB[e])-mean(pA[e]))-(mean(pB[!e])-mean(pA[!e]))}

specs <- list(A=c("g13_5"), MEAN=c("gmean_5"), MEAN_PLUS_G13=c("g13_5","gmean_5"),
              B1=c("g13_5","low11"), B2=c("g13_5","grp4"), B3=c("g13_5","g11_5"))
for (outcome in c(1,2)) {
  lab <- c("disability","death")[outcome]
  dd <- d0[d0$oc %in% c(0,outcome), ]; y <- as.numeric(dd$oc==outcome)
  cat(sprintf("\n--- %s: n=%d events=%d ---\n", lab, nrow(dd), sum(y)))
  fits <- list(); preds <- list()
  for (nm in names(specs)) {
    f <- as.formula(paste("y ~", paste(c(COV, specs[[nm]]), collapse="+")))
    fits[[nm]] <- glm(f, data=dd, family=binomial()); preds[[nm]] <- predict(fits[[nm]], type="response")
  }
  rA <- roc(y, preds$A, quiet=TRUE)
  cf <- summary(fits$A)$coefficients
  cat(sprintf("  OR per 5 kg (model A, 2013 grip): %.3f (%.3f-%.3f) p=%.3g\n",
              exp(cf["g13_5","Estimate"]), exp(cf["g13_5","Estimate"]-1.96*cf["g13_5","Std. Error"]),
              exp(cf["g13_5","Estimate"]+1.96*cf["g13_5","Std. Error"]), cf["g13_5","Pr(>|z|)"]))
  cat(sprintf("  AUC A=%.4f | MEAN=%.4f | MEAN+G13=%.4f | B2=%.4f\n",
              as.numeric(auc(rA)), as.numeric(auc(roc(y,preds$MEAN,quiet=TRUE))),
              as.numeric(auc(roc(y,preds$MEAN_PLUS_G13,quiet=TRUE))),
              as.numeric(auc(roc(y,preds$B2,quiet=TRUE)))))
  for (nm in c("MEAN","MEAN_PLUS_G13","B2")) {
    dl <- roc.test(rA, roc(y,preds[[nm]],quiet=TRUE), method="delong")
    lr <- anova(fits$A, fits[[nm]], test="Chisq")
    cat(sprintf("  %-14s dAUC=%+.4f DeLong p=%.3f dAIC=%+.2f LR chi2=%.2f df=%d p=%.3f cNRI=%+.3f IDI=%.4f\n",
                nm, as.numeric(auc(roc(y,preds[[nm]],quiet=TRUE)))-as.numeric(auc(rA)),
                dl$p.value, AIC(fits[[nm]])-AIC(fits$A), lr$Deviance[2], lr$Df[2], lr$`Pr(>Chi)`[2],
                nri_cont(preds$A,preds[[nm]],y), idi(preds$A,preds[[nm]],y)))
  }
  # cluster bootstrap CI for dAUC of MEAN and MEAN_PLUS_G13
  cl <- unique(dd$communityID); B <- 300
  store <- matrix(NA, B, 2, dimnames=list(NULL, c("MEAN","MEAN_PLUS_G13")))
  for (i in 1:B) {
    pick <- sample(cl, length(cl), replace=TRUE)
    idx <- unlist(lapply(pick, function(k) which(dd$communityID==k)))
    d2 <- dd[idx,]; y2 <- y[idx]
    if (length(unique(y2))<2) next
    pA2 <- predict(glm(as.formula(paste("y ~",paste(c(COV,"g13_5"),collapse="+"))),data=d2,family=binomial()),type="response")
    aA2 <- as.numeric(auc(roc(y2,pA2,quiet=TRUE)))
    for (nm in c("MEAN","MEAN_PLUS_G13")) {
      pX <- tryCatch(predict(glm(as.formula(paste("y ~",paste(c(COV,specs[[nm]]),collapse="+"))),data=d2,family=binomial()),type="response"), error=function(e) NULL)
      if (!is.null(pX)) store[i,nm] <- as.numeric(auc(roc(y2,pX,quiet=TRUE))) - aA2
    }
  }
  for (nm in c("MEAN","MEAN_PLUS_G13")) {
    v <- store[,nm]; v <- v[is.finite(v)]
    cat(sprintf("  bootstrap dAUC %s: 95%%CI %+.4f to %+.4f (B=%d)\n", nm, quantile(v,.025), quantile(v,.975), length(v)))
  }
}

# ---------------------------------------------------------------
cat("\n=============================================================\n")
cat("(D) included vs excluded for missing 2013 ADL (review M5)\n")
cat("=============================================================\n")
two <- m %>% filter(!is.na(g11), !is.na(g13), !is.na(grp19), !is.na(age13), !is.na(male))
cat("two-wave classifiable n =", nrow(two), "\n")
cat("a13 distribution:\n"); print(table(two$a13, useNA="ifany"))
two$excluded <- as.integer(!(two$a13 %in% 0))
cmp <- function(v, lab, digits=1) {
  a <- mean(two[[v]][two$excluded==0], na.rm=TRUE); b <- mean(two[[v]][two$excluded==1], na.rm=TRUE)
  cat(sprintf("  %-28s included %6.*f | excluded %6.*f\n", lab, digits, a, digits, b))
}
cmp("age13","age, years"); cmp("male","male, %",0); 
two$malep <- 100*two$male
cmp("malep","male, %",1)
cmp("g13","grip 2013, kg"); cmp("g11","grip 2011, kg")
cmp("bmi","BMI"); cmp("nchronic","chronic conditions"); cmp("cesd","CES-D-10")
cmp("low11_19","low grip 2011, %",0)
two$low11p <- 100*two$low11_19; cmp("low11p","low grip 2011, %",1)
cat("  missing-data counts among the excluded: bmi", sum(is.na(two$bmi[two$excluded==1])),
    "| cesd", sum(is.na(two$cesd[two$excluded==1])), "\n")
cat("  death by 2018: included", sprintf("%.1f%%", 100*mean(two$dead18[two$excluded==0], na.rm=TRUE)),
    "| excluded", sprintf("%.1f%%", 100*mean(two$dead18[two$excluded==1], na.rm=TRUE)), "\n")

sink()
cat("\nDONE\n")
