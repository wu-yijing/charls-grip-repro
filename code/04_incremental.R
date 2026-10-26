# 04_incremental.R -- Analysis (3): incremental predictive value of the 2011 reading (reviewer M6)
# Model A = covariates + 2013 grip (current reading only)
# Model B = Model A + 2011 status (pattern / prior low-grip flag)
# Reports LR test, delta-AIC, AUC (DeLong), continuous NRI (bootstrap CI), IDI.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
m <- readRDS("out/jamda_full.rds")
set.seed(20261023)
sink("out/S3_incremental_report.txt", split=TRUE)
cat("==================================================================\n")
cat("SUPPLEMENTARY TABLE S3 (generated", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
cat("Analysis 3 (M6): does the 2011 reading add anything to the 2013 reading?\n")
cat("==================================================================\n\n")

b <- get_base(m, "grp19")
b$g13_5  <- b$g13/5
b$g11_5  <- b$g11/5
b$low13  <- b$low13_19
b$low11  <- b$low11_19
b$grp4   <- b$grp
d0 <- b[complete.cases(b[, c(COV, "g13_5", "g11_5", "low11", "low13", "grp4", "oc")]), ]
cat("analysis set (complete covariates + both grips): n =", nrow(d0), "\n")
cat("outcomes: alive-free", sum(d0$oc==0), "| disabled", sum(d0$oc==1), "| dead", sum(d0$oc==2), "\n\n")

nri_cont <- function(pA, pB, event){
  up <- pB > pA; dn <- pB < pA
  e <- event==1
  nri_e  <- mean(up[e])  - mean(dn[e])
  nri_ne <- mean(dn[!e]) - mean(up[!e])
  c(NRI=nri_e+nri_ne, NRI_event=nri_e, NRI_nonevent=nri_ne,
    up_event=mean(up[e]), down_event=mean(dn[e]),
    up_nonevent=mean(up[!e]), down_nonevent=mean(dn[!e]))
}
idi <- function(pA, pB, event){ e <- event==1
  mean(pB[e])-mean(pA[e]) - (mean(pB[!e])-mean(pA[!e])) }

run <- function(outcome, label){
  dd <- d0[d0$oc %in% c(0, outcome), ]
  y <- as.numeric(dd$oc==outcome)
  fA  <- as.formula(paste("y ~", paste(c(COV,"g13_5"), collapse="+")))
  fB1 <- as.formula(paste("y ~", paste(c(COV,"g13_5","low11"), collapse="+")))
  fB2 <- as.formula(paste("y ~", paste(c(COV,"g13_5","grp4"), collapse="+")))
  fA0 <- as.formula(paste("y ~", paste(c(COV,"g11_5"), collapse="+")))  # 2011 only, for comparison
  gA  <- glm(fA,  data=dd, family=binomial());  gB1 <- glm(fB1, data=dd, family=binomial())
  gB2 <- glm(fB2, data=dd, family=binomial()); gA0 <- glm(fA0, data=dd, family=binomial())
  pA <- predict(gA, type="response"); pB1 <- predict(gB1, type="response"); pB2 <- predict(gB2, type="response")
  rA <- roc(y, pA, quiet=TRUE); rB1 <- roc(y, pB1, quiet=TRUE); rB2 <- roc(y, pB2, quiet=TRUE)
  d1 <- roc.test(rA, rB1, method="delong"); d2 <- roc.test(rA, rB2, method="delong")
  lr1 <- anova(gA, gB1, test="Chisq"); lr2 <- anova(gA, gB2, test="Chisq")
  n1 <- nri_cont(pA, pB1, y); n2 <- nri_cont(pA, pB2, y)
  cat("------ outcome:", label, "(n =", nrow(dd), ", events =", sum(y), ") ------\n")
  cat(sprintf("  AUC 2013 grip only (A)   : %.3f (%.3f-%.3f)\n", as.numeric(auc(rA)),
      as.numeric(ci.auc(rA))[1], as.numeric(ci.auc(rA))[3]))
  cat(sprintf("  AUC A + 2011 low flag (B1): %.3f (%.3f-%.3f) | DeLong p = %.3f | dAUC = %+.4f\n",
      as.numeric(auc(rB1)), as.numeric(ci.auc(rB1))[1], as.numeric(ci.auc(rB1))[3],
      d1$p.value, as.numeric(auc(rB1))-as.numeric(auc(rA))))
  cat(sprintf("  AUC A + 4-group pattern(B2): %.3f (%.3f-%.3f) | DeLong p = %.3f | dAUC = %+.4f\n",
      as.numeric(auc(rB2)), as.numeric(ci.auc(rB2))[1], as.numeric(ci.auc(rB2))[3],
      d2$p.value, as.numeric(auc(rB2))-as.numeric(auc(rA))))
  cat(sprintf("  AUC 2011 grip only       : %.3f\n", as.numeric(auc(roc(y, predict(gA0, type='response'), quiet=TRUE)))))
  cat(sprintf("  LR test  B1 vs A: chi2(%d) = %.2f, p = %.4f ; dAIC = %+.2f\n",
      lr1$Df[2], lr1$Deviance[2], lr1$`Pr(>Chi)`[2], AIC(gB1)-AIC(gA)))
  cat(sprintf("  LR test  B2 vs A: chi2(%d) = %.2f, p = %.4f ; dAIC = %+.2f\n",
      lr2$Df[2], lr2$Deviance[2], lr2$`Pr(>Chi)`[2], AIC(gB2)-AIC(gA)))
  cat(sprintf("  continuous NRI  B1 vs A: %+.3f (events %+.3f / non-events %+.3f)\n", n1[1], n1[2], n1[3]))
  cat(sprintf("  continuous NRI  B2 vs A: %+.3f (events %+.3f / non-events %+.3f)\n", n2[1], n2[2], n2[3]))
  cat(sprintf("  IDI B2 vs A: %+.4f | B1 vs A: %+.4f\n", idi(pA,pB2,y), idi(pA,pB1,y)))
  # cluster bootstrap for NRI (B2) and dAUC
  B <- 500; nb <- length(y); ids <- dd$communityID
  cl <- unique(ids); nri_b <- numeric(B); dauc_b <- numeric(B)
  for (i in 1:B) {
    pick <- sample(cl, length(cl), replace=TRUE)
    idx <- unlist(lapply(pick, function(k) which(ids==k)))
    dd2 <- dd[idx, ]; y2 <- y[idx]
    if (length(unique(y2))<2) { nri_b[i] <- NA; dauc_b[i] <- NA; next }
    gA2 <- glm(fA, data=dd2, family=binomial()); gB2_ <- glm(fB2, data=dd2, family=binomial())
    pA2 <- predict(gA2, type="response"); pB2_ <- predict(gB2_, type="response")
    nri_b[i] <- nri_cont(pA2, pB2_, y2)[1]
    dauc_b[i] <- as.numeric(auc(roc(y2, pB2_, quiet=TRUE)) - auc(roc(y2, pA2, quiet=TRUE)))
  }
  cat(sprintf("  bootstrap (500 cluster reps): NRI(B2) 95%%CI %.3f to %.3f | dAUC 95%%CI %+.4f to %+.4f\n\n",
      quantile(nri_b, .025, na.rm=TRUE), quantile(nri_b, .975, na.rm=TRUE),
      quantile(dauc_b, .025, na.rm=TRUE), quantile(dauc_b, .975, na.rm=TRUE)))
  data.frame(outcome=label, n=nrow(dd), events=sum(y),
    auc_A=as.numeric(auc(rA)), auc_B1=as.numeric(auc(rB1)), auc_B2=as.numeric(auc(rB2)),
    p_delong_B1=d1$p.value, p_delong_B2=d2$p.value,
    dAIC_B1=AIC(gB1)-AIC(gA), dAIC_B2=AIC(gB2)-AIC(gA),
    LR_B1=lr1$Deviance[2], p_LR_B1=lr1$`Pr(>Chi)`[2],
    LR_B2=lr2$Deviance[2], p_LR_B2=lr2$`Pr(>Chi)`[2],
    NRI_B1=n1[1], NRI_B2=n2[1], NRI_B2_lo=quantile(nri_b,.025,na.rm=TRUE),
    NRI_B2_hi=quantile(nri_b,.975,na.rm=TRUE), IDI_B2=idi(pA,pB2,y))
}
res <- rbind(run(1, "disability (2018)"), run(2, "death (2018)"))
write.csv(res, "out/S3_incremental.csv", row.names=FALSE)
cat("--- summary table ---\n"); print(res, row.names=FALSE, digits=3)

## ---- sensitivity: AWGS 2025 cut-offs for the 4-group pattern ----
cat("\n--- sensitivity: four-group pattern defined by AWGS 2025 cut-offs ---\n")
b25 <- get_base(m, "grp25")
b25$g13_5 <- b25$g13/5; b25$low11 <- b25$low11_25; b25$grp4 <- b25$grp
d25 <- b25[complete.cases(b25[, c(COV,"g13_5","low11","grp4","oc")]), ]
for (outcome in c(1,2)) {
  dd <- d25[d25$oc %in% c(0,outcome), ]; y <- as.numeric(dd$oc==outcome)
  fA <- as.formula(paste("y ~", paste(c(COV,"g13_5"), collapse="+")))
  fB <- as.formula(paste("y ~", paste(c(COV,"g13_5","grp4"), collapse="+")))
  gA <- glm(fA, data=dd, family=binomial()); gB <- glm(fB, data=dd, family=binomial())
  rA <- roc(y, predict(gA,type="response"), quiet=TRUE); rB <- roc(y, predict(gB,type="response"), quiet=TRUE)
  lr <- anova(gA, gB, test="Chisq")
  cat(sprintf("  %-18s n=%d | AUC %.3f -> %.3f | DeLong p=%.3f | dAIC=%+.2f | LR p=%.3f\n",
    c("disability","death")[outcome], nrow(dd), as.numeric(auc(rA)), as.numeric(auc(rB)),
    roc.test(rA,rB,method="delong")$p.value, AIC(gB)-AIC(gA), lr$`Pr(>Chi)`[2]))
}
sink()
cat("\nS3 done -> out/S3_incremental.csv, out/S3_incremental_report.txt\n")
