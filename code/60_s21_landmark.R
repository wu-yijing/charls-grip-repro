# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# Supplementary Table S21: incremental value of the 2011 reading once the 2013
# ADL-free landmark restriction is removed.
# 60_s21_landmark.R -- S21: incremental value of the 2011 reading when the
# 2013 "ADL-free" landmark restriction is removed (reviewer M1).
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")

m <- readRDS("out/jamda_full.rds")
set.seed(20261028)

b <- m %>% filter(!is.na(grp19), !is.na(age13), !is.na(male))
b$grp <- relevel(factor(b$grp19, levels=GRP), ref="never")
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5; b$gmean_5 <- (b$g11 + b$g13)/10
b$low11 <- b$low11_19; b$grp4 <- b$grp
need <- c(COV, "g13_5","g11_5","gmean_5","low11","grp4","oc","a13")
b <- b[complete.cases(b[, need]), ]

specs <- list(A = "g13_5", B1 = c("g13_5","low11"),
              B2 = c("g13_5","grp4"), B3 = c("g13_5","g11_5"))
fit1 <- function(d, vy, xs)
  glm(as.formula(paste(vy, "~", paste(c(COV, xs), collapse="+"))), data=d, family=binomial())

rows <- list(); orows <- list()
analyse <- function(d, vy, lab){
  y <- as.numeric(d[[vy]])
  fits <- lapply(specs, function(s) fit1(d, vy, s))
  pr <- lapply(fits, function(f) predict(f, type="response"))
  rA <- roc(y, pr$A, quiet=TRUE)
  for (nm in names(specs)[-1]) {
    rX <- roc(y, pr[[nm]], quiet=TRUE)
    lf <- anova(fits$A, fits[[nm]], test="Chisq")
    rows[[length(rows)+1]] <<- data.frame(set=lab, model=nm, n=nrow(d), events=sum(y),
      auc_A=round(as.numeric(auc(rA)),4), auc_X=round(as.numeric(auc(rX)),4),
      dAUC=round(as.numeric(auc(rX))-as.numeric(auc(rA)),5),
      p_Delong=signif(roc.test(rA, rX, method="delong")$p.value,3),
      LR=round(lf$Deviance[2],2), df=lf$Df[2], p_LR=signif(lf$`Pr(>Chi)`[2],3),
      dAIC=round(AIC(fits[[nm]])-AIC(fits$A),2))
  }
  f6 <- fit1(d, vy, "gmean_5"); p6 <- predict(f6, type="response"); r6 <- roc(y, p6, quiet=TRUE)
  rows[[length(rows)+1]] <<- data.frame(set=lab, model="B6", n=nrow(d), events=sum(y),
    auc_A=round(as.numeric(auc(rA)),4), auc_X=round(as.numeric(auc(r6)),4),
    dAUC=round(as.numeric(auc(r6))-as.numeric(auc(rA)),5),
    p_Delong=signif(roc.test(rA, r6, method="delong")$p.value,3),
    LR=NA, df=NA, p_LR=NA, dAIC=round(AIC(f6)-AIC(fits$A),2))
  cf <- summary(fits$B2)$coefficients
  for (g in c("remitted","incident","persistent")) {
    k <- paste0("grp4", g)
    if (k %in% rownames(cf))
      orows[[length(orows)+1]] <<- data.frame(set=lab, group=g, n=nrow(d), events=sum(y),
        OR=round(exp(cf[k,1]),2), lo=round(exp(cf[k,1]-1.96*cf[k,2]),2),
        hi=round(exp(cf[k,1]+1.96*cf[k,2]),2), p=signif(cf[k,4],3))
  }
}

for (samp in c("restricted","unrestricted")) {
  d0 <- if (samp=="restricted") b[b$a13 == 0, ] else b[b$a13 %in% c(0,1), ]
  d0$y_dis <- as.numeric(d0$oc == 1)
  d0$y_dea <- as.numeric(d0$oc == 2)
  d0$y_cmp <- as.numeric(d0$oc %in% c(1,2))
  analyse(d0[d0$oc %in% c(0,1), ], "y_dis", paste0(samp, " | disability"))
  analyse(d0[d0$oc %in% c(0,2), ], "y_dea", paste0(samp, " | death"))
  analyse(d0, "y_cmp", paste0(samp, " | disability or death"))
}
res <- do.call(rbind, rows); ors <- do.call(rbind, orows)
print(res, row.names=FALSE)
cat("\n=== four-group ORs (B2, reference never low) ===\n")
print(ors, row.names=FALSE)
write.csv(res, "out/S21_landmark.csv", row.names=FALSE)
write.csv(ors, "out/S21_landmark_OR.csv", row.names=FALSE)
