# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# Supplementary Table S21 Panel C: unrestricted sample with 2013 ADL status in each model.
# 63_s21_panelC.R -- S21 Panel C: unrestricted sample, 2013 ADL status in every model
suppressMessages(library(pROC)); 
m <- readRDS("out/jamda_full.rds")
COV <- c("age10","male","bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4")
b <- m[!is.na(m$grp19) & !is.na(m$age13) & !is.na(m$male), ]
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5; b$gmean_5 <- (b$g11+b$g13)/10
b$low11 <- b$low11_19; b$grp4 <- factor(b$grp19, levels=c("never","remitted","incident","persistent"))
b <- b[complete.cases(b[, c(COV,"g13_5","g11_5","gmean_5","low11","grp4","oc","a13")]), ]
b <- b[b$a13 %in% c(0,1), ]
b$y_dis <- as.numeric(b$oc==1); b$y_dea <- as.numeric(b$oc==2); b$y_cmp <- as.numeric(b$oc %in% c(1,2))

A <- c(COV, "a13")
fit <- function(d,v,xs) glm(as.formula(paste(v,"~",paste(c(A,xs),collapse="+"))), data=d, family=binomial())
specs <- list(B1=c("g13_5","low11"), B2=c("g13_5","grp4"), B3=c("g13_5","g11_5"))
rows <- list()
for (v in c("y_dis","y_dea","y_cmp")) {
  d <- switch(v, y_dis=b[b$oc %in% c(0,1),], y_dea=b[b$oc %in% c(0,2),], y_cmp=b)
  fA <- fit(d,v,c("g13_5")); pA <- predict(fA,type="response"); rA <- roc(d[[v]],pA,quiet=TRUE)
  for (nm in names(specs)) {
    fX <- fit(d,v,specs[[nm]]); pX <- predict(fX,type="response"); rX <- roc(d[[v]],pX,quiet=TRUE)
    lf <- anova(fA,fX,test="Chisq")
    rows[[length(rows)+1]] <- data.frame(outcome=v, model=nm, n=nrow(d), ev=sum(d[[v]]),
      aucA=round(as.numeric(auc(rA)),3), aucX=round(as.numeric(auc(rX)),3),
      dAUC=round(as.numeric(auc(rX))-as.numeric(auc(rA)),5),
      P=signif(roc.test(rA,rX,method="delong")$p.value,3),
      LR=round(lf$Deviance[2],2), df=lf$Df[2], pLR=signif(lf$`Pr(>Chi)`[2],3),
      dAIC=round(AIC(fX)-AIC(fA),2))
  }
  # B6: mean of two readings replaces current value (not nested)
  f6 <- fit(d,v,"gmean_5"); p6 <- predict(f6,type="response"); r6 <- roc(d[[v]],p6,quiet=TRUE)
  rows[[length(rows)+1]] <- data.frame(outcome=v, model="B6", n=nrow(d), ev=sum(d[[v]]),
    aucA=round(as.numeric(auc(rA)),3), aucX=round(as.numeric(auc(r6)),3),
    dAUC=round(as.numeric(auc(r6))-as.numeric(auc(rA)),5),
    P=signif(roc.test(rA,r6,method="delong")$p.value,3), LR=NA, df=NA, pLR=NA,
    dAIC=round(AIC(f6)-AIC(fA),2))
}
r <- do.call(rbind, rows)
for (i in seq_len(nrow(r))) cat(sprintf("%-8s %-2s n=%4d ev=%3d AUC=%.3f dAUC=%+.5f P=%.3f LR=%.2f df=%s pLR=%s dAIC=%+.2f\n",
  r$outcome[i], r$model[i], r$n[i], r$ev[i], r$aucA[i], r$dAUC[i], r$P[i], r$LR[i],
  as.character(r$df[i]), as.character(r$pLR[i]), r$dAIC[i]))
write.csv(r, "out/S21_panelC.csv", row.names=FALSE)
