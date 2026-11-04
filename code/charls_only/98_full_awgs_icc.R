# 98_full_awgs_icc.R -- Supplementary Tables S37 and S38 (archived with v1.0.8)
# Run from the analysis directory that contains out/jamda_full.rds, _x/ and lib.R;
# see code/charls_only/README.md. No machine-specific paths are used.
# (1) full AWGS 2019 "possible sarcopenia" definition in the whole 6,052 analytic sample
# (2) short-interval ICC 0.80-0.99 sensitivity table
# All inputs are the ORIGINAL local CHARLS files; no data are simulated.
suppressMessages({library(haven); library(dplyr)})
source("lib.R")
out_dir <- "out"

m <- readRDS("out/jamda_full.rds") %>% mutate(ID = as.character(ID))
# analytic sample, identical to jfa_00_persist_boot.R / jfa_08_external_anchor.R
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
cat(sprintf("analytic n = %d | communities = %d\n", nrow(b), length(unique(b$communityID))))

## ---------------- (1) full AWGS 2019 possible sarcopenia ----------------
fix11 <- function(x) paste0(substr(x, 1, 10), "0", substr(x, 11, 11))
a <- read_dta("_x/2011/biomarkers.dta") %>% mutate(ID = fix11(as.character(ID))) %>%
  transmute(ID, cs1 = suppressWarnings(as.numeric(qh003)))
c13 <- read_dta("_x/2013/Biomarker.dta") %>% mutate(ID = as.character(ID)) %>%
  transmute(ID, cs1b = suppressWarnings(as.numeric(qh003)))

d <- b %>% select(ID, g11, g13, age11, age13, male, communityID, low11_19, low13_19) %>%
  left_join(a, by = "ID") %>% left_join(c13, by = "ID") %>%
  mutate(cs11 = ifelse(!is.na(cs1)  & cs1  > 0 & cs1  < 120, cs1,  NA),
         cs13 = ifelse(!is.na(cs1b) & cs1b > 0 & cs1b < 120, cs1b, NA),
         slow11 = !is.na(cs11) & cs11 >= 12,
         slow13 = !is.na(cs13) & cs13 >= 12,
         grip11 = low11_19 == 1,
         grip13 = low13_19 == 1,
         poss11 = grip11 | slow11,
         poss13 = grip13 | slow13)

cat(sprintf("chair stand available 2011: %d | 2013: %d | both: %d | neither-pair missing: %d\n",
            sum(!is.na(d$cs11)), sum(!is.na(d$cs13)),
            sum(!is.na(d$cs11) & !is.na(d$cs13)),
            sum(is.na(d$cs11) | is.na(d$cs13))))

## cluster bootstrap for a 2x2 table (same machinery as jfa_00_persist_boot.R)
prep2 <- function(dd, x11, x13){
  cl <- as.integer(factor(dd$communityID))
  A <- tapply(as.integer(dd[[x11]] & dd[[x13]]), cl, sum, default=0)
  Bv<- tapply(as.integer(dd[[x11]] & !dd[[x13]]), cl, sum, default=0)
  C <- tapply(as.integer(!dd[[x11]] & dd[[x13]]), cl, sum, default=0)
  D <- tapply(as.integer(!dd[[x11]] & !dd[[x13]]), cl, sum, default=0)
  list(A=as.numeric(A),B=as.numeric(Bv),C=as.numeric(C),D=as.numeric(D),
       G=length(unique(cl)), cells=c(sum(A),sum(Bv),sum(C),sum(D)))
}
idx2 <- function(A,Bv,C,D){
  N <- A+Bv+C+D; po <- (A+D)/N
  pe <- ((A+Bv)*(A+C) + (C+D)*(Bv+D))/N^2
  c(persist=A/(A+Bv), revert=Bv/(A+Bv), incid=C/(C+D),
    pa=2*A/(2*A+Bv+C), kappa=(po-pe)/(1-pe), pabak=2*po-1,
    prev=(A+Bv)/N, n=N)
}
bootci <- function(p, B=2000, seed=20261104){
  set.seed(seed); o <- matrix(NA_real_, B, 8)
  gi <- seq_len(p$G)
  for(i in 1:B){
    w <- tabulate(sample(gi, p$G, replace=TRUE), nbins=p$G)
    o[i,] <- as.numeric(idx2(sum(w*p$A), sum(w*p$B), sum(w*p$C), sum(w*p$D)))
  }
  lo <- apply(o, 2, quantile, 0.025, na.rm=TRUE, names=FALSE)
  hi <- apply(o, 2, quantile, 0.975, na.rm=TRUE, names=FALSE)
  rbind(lo=lo, hi=hi)
}
report2 <- function(dd, x11, x13, tag){
  p <- prep2(dd, x11, x13)
  est <- idx2(p$cells[1],p$cells[2],p$cells[3],p$cells[4])
  ci  <- bootci(p)
  data.frame(definition=tag, metric=names(est), est=as.numeric(est),
             lo=ci[1,], hi=ci[2,], n=sum(p$cells), clusters=p$G,
             cells=paste(p$cells, collapse="/"), row.names=NULL)
}
o_grip <- report2(d, "grip11", "grip13", "Low grip strength only (primary)")
o_full <- report2(d, "poss11", "poss13", "Full AWGS 2019 possible sarcopenia (low grip OR chair stand >=12 s)")
o_both <- rbind(o_grip, o_full)
cat("\n=== full-sample comparison ===\n")
print(o_both[o_both$metric %in% c("persist","revert","pa","kappa","prev"),], row.names=FALSE, digits=4)
write.csv(o_both, file.path(out_dir,"JFA_full_awgs.csv"), row.names=FALSE)

cat("\n=== 2011 classification counts (full definition) ===\n")
p2 <- d %>% mutate(type=case_when(grip11 & slow11 ~ "both", grip11 ~ "grip only", TRUE ~ "chair only"))
tab <- p2 %>% filter(poss11) %>% group_by(type) %>%
  summarise(n=n(), not_reproduced=sum(!poss13), pct=100*mean(!poss13), .groups="drop")
print(as.data.frame(tab))
write.csv(tab, file.path(out_dir,"JFA_full_awgs_decomp.csv"), row.names=FALSE)

## ---------------- (2) ICC 0.80 - 0.99 sensitivity ----------------
nv <- b[b$grp19=="never",]
Phi2 <- function(h,k,r) integrate(function(x) dnorm(x)*pnorm((k-r*x)/sqrt(1-r^2)), -Inf, h, rel.tol=1e-9)$value
exp_rev <- function(p,r){ cc<-qnorm(p); 1 - Phi2(cc,cc,r)/p }
sd_m <- sd(b$g11[b$male==1]); sd_f <- sd(b$g11[b$male==0])
p_m  <- mean(b$low11_19[b$male==1]); p_f <- mean(b$low11_19[b$male==0])
n_m  <- sum(b$low11_19==1 & b$male==1); n_f <- sum(b$low11_19==1 & b$male==0)
sd_chg <- sd(nv$g13-nv$g11); var_d <- sd_chg^2
cat(sprintf("\nSD(change)=%.3f | SD grip men=%.3f women=%.3f | p low men=%.4f women=%.4f\n",
            sd_chg, sd_m, sd_f, p_m, p_f))

rows <- list()
for (icc in c(0.80,0.85,0.90,0.95,0.99)){
  rev_m <- exp_rev(p_m,icc); rev_f <- exp_rev(p_f,icc)
  rev <- (rev_m*n_m + rev_f*n_f)/(n_m+n_f)
  se_m <- sd_m*sqrt(1-icc); se_f <- sd_f*sqrt(1-icc)
  se_pool <- sqrt(mean(c(se_m^2,se_f^2)))
  var_true <- max(var_d - 2*se_pool^2, 0)
  rows[[length(rows)+1]] <- data.frame(
    icc=icc, error_SD=se_pool, expected_reversion=100*rev,
    obs_over_exp=59.505/(100*rev),
    var_share_error=100*2*se_pool^2/var_d,
    true_change_SD=sqrt(var_true))
}
sens <- do.call(rbind, rows)
cat("\n=== ICC 0.80-0.99 sensitivity ===\n"); print(sens, row.names=FALSE, digits=4)
write.csv(sens, file.path(out_dir,"JFA_icc_sensitivity.csv"), row.names=FALSE)
cat("\nwritten: JFA_full_awgs.csv, JFA_full_awgs_decomp.csv, JFA_icc_sensitivity.csv\n")
