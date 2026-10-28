# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# Supplementary Table S21 Panel C: 2013 ADL status added to every unrestricted model.
# 62_s21_adjust.R -- does the S21 increase survive adjustment for 2013 ADL status?
suppressMessages(library(pROC)); 
m <- readRDS("out/jamda_full.rds")
COV <- c("age10","male","bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4")
b <- m[!is.na(m$grp19) & !is.na(m$age13) & !is.na(m$male), ]
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5
b$low11 <- b$low11_19; b$grp4 <- factor(b$grp19, levels=c("never","remitted","incident","persistent"))
b <- b[complete.cases(b[, c(COV,"g13_5","g11_5","low11","grp4","oc","a13")]), ]
b <- b[b$a13 %in% c(0,1), ]
b$y_dis <- as.numeric(b$oc == 1); b$y_cmp <- as.numeric(b$oc %in% c(1,2)); b$y_dea <- as.numeric(b$oc == 2)

fit <- function(d, y, xs) glm(as.formula(paste(y,"~",paste(xs,collapse="+"))), data=d, family=binomial())
run <- function(d, v, extra, lab){
  d$y <- as.numeric(d[[v]])
  A <- COV; if (extra) A <- c(COV, "a13")
  fA <- fit(d,"y",c(A,"g13_5")); pA <- predict(fA,type="response")
  out <- c()
  for (nm in c("B2","B3")) {
    xs <- if (nm=="B2") c(A,"g13_5","grp4") else c(A,"g13_5","g11_5")
    fX <- fit(d,"y",xs); pX <- predict(fX,type="response")
    rA <- roc(d$y,pA,quiet=TRUE); rX <- roc(d$y,pX,quiet=TRUE)
    out <- c(out, sprintf("%s dAUC=%+.4f P=%.3f", nm,
      as.numeric(auc(rX))-as.numeric(auc(rA)), roc.test(rA,rX,method="delong")$p.value))
  }
  cat(sprintf("%-42s n=%d ev=%d | %s\n", lab, nrow(d), sum(d$y), paste(out, collapse=" | ")))
}
for (extra in c(FALSE, TRUE)) {
  tag <- if (extra) "+2013 ADL as covariate" else "no adjustment for 2013 ADL"
  run(b[b$oc %in% c(0,1),], "y_dis", extra, paste("unrestricted disability", tag))
  run(b, "y_cmp", extra, paste("unrestricted composite ", tag))
  run(b[b$oc %in% c(0,2),], "y_dea", extra, paste("unrestricted death     ", tag))
}
