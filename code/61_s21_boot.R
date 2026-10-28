# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds. Nothing here reads CHARLS data.
#
# Supplementary Table S21: cluster bootstrap for the change in AUC in the unrestricted sample.
# 61_s21_boot.R -- cluster bootstrap for S21 dAUC (unrestricted sample)
suppressMessages({library(pROC)})

m <- readRDS("out/jamda_full.rds"); set.seed(20261028)
COV <- c("age10","male","bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4")
b <- m[!is.na(m$grp19) & !is.na(m$age13) & !is.na(m$male), ]
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5; b$gmean_5 <- (b$g11+b$g13)/10
b$low11 <- b$low11_19; b$grp4 <- factor(b$grp19, levels=c("never","remitted","incident","persistent"))
b <- b[complete.cases(b[, c(COV,"g13_5","g11_5","gmean_5","low11","grp4","oc","a13","communityID")]), ]
b <- b[b$a13 %in% c(0,1), ]
b$y_dis <- as.numeric(b$oc == 1); b$y_dea <- as.numeric(b$oc == 2); b$y_cmp <- as.numeric(b$oc %in% c(1,2))

auc_of <- function(d, y, xs){
  f <- glm(as.formula(paste("y ~", paste(c(COV, xs), collapse="+"))), data=d, family=binomial())
  as.numeric(auc(roc(y, predict(f, type="response"), quiet=TRUE)))
}
run <- function(d, v, specs, B=200){
  d$y <- as.numeric(d[[v]])
  idx_by_cl <- split(seq_len(nrow(d)), d$communityID)
  cl <- names(idx_by_cl)
  base <- auc_of(d, d$y, "g13_5")
  res <- sapply(names(specs), function(nm){
    v2 <- numeric(B)
    for (i in seq_len(B)) {
      pick <- sample(cl, length(cl), TRUE)
      ii <- unlist(idx_by_cl[pick], use.names=FALSE)
      d2 <- d[ii, ]
      if (length(unique(d2$y)) < 2) { v2[i] <- NA; next }
      v2[i] <- tryCatch(auc_of(d2, d2$y, specs[[nm]]) - auc_of(d2, d2$y, "g13_5"), error=function(e) NA)
    }
    v2 <- v2[is.finite(v2)]
    sprintf("%+.4f (%+.4f to %+.4f) B=%d", mean(v2), quantile(v2,.025), quantile(v2,.975), length(v2))
  })
  cat(sprintf("%-22s baseAUC=%.3f | %s\n", v, base, paste(names(res), res, sep="=", collapse=" | ")))
}
specs <- list(B2=c("g13_5","grp4"), B3=c("g13_5","g11_5"))
run(b[b$oc %in% c(0,1), ], "y_dis", specs)
run(b, "y_cmp", specs)
run(b[b$oc %in% c(0,2), ], "y_dea", specs)
