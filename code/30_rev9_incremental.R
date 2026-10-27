# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 30_rev9_incremental.R -- rev9 addition for reviewer M4-1
# Purpose: the published incremental analysis only entered the prior reading as a
# BINARY flag / 4-group pattern. Add CONTINUOUS prior reading and change score.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
m <- readRDS("out/jamda_full.rds")
sink("out/30_rev9_incremental.txt", split = TRUE)
set.seed(20261024)

b <- get_base(m, "grp19")
b$g13_5 <- b$g13/5
b$g11_5 <- b$g11/5
b$chg_5 <- (b$g13 - b$g11)/5
b$low11 <- b$low11_19
b$low13 <- b$low13_19
b$grp4  <- b$grp
d0 <- b[complete.cases(b[, c(COV,"g13_5","g11_5","chg_5","low11","low13","grp4","oc")]), ]
cat("analysis set n =", nrow(d0), "\n")

nri_cont <- function(pA,pB,event){up<-pB>pA;dn<-pB<pA;e<-event==1
  (mean(up[e])-mean(dn[e]))+(mean(dn[!e])-mean(up[!e]))}
idi <- function(pA,pB,event){e<-event==1
  (mean(pB[e])-mean(pA[e]))-(mean(pB[!e])-mean(pA[!e]))}

specs <- list(
  A   = c("g13_5"),
  B1  = c("g13_5","low11"),
  B2  = c("g13_5","grp4"),
  B3  = c("g13_5","g11_5"),
  B4  = c("g13_5","chg_5"),
  B5  = c("g13_5","g11_5","low11")
)

out <- list()
for (outcome in c(1,2)) {
  lab <- c("disability","death")[outcome]
  dd <- d0[d0$oc %in% c(0,outcome), ]
  y  <- as.numeric(dd$oc==outcome)
  fits <- list(); preds <- list()
  for (nm in names(specs)) {
    f <- as.formula(paste("y ~", paste(c(COV, specs[[nm]]), collapse="+")))
    fits[[nm]] <- glm(f, data=dd, family=binomial())
    preds[[nm]] <- predict(fits[[nm]], type="response")
  }
  rA <- roc(y, preds$A, quiet=TRUE)
  for (nm in names(specs)[-1]) {
    rX <- roc(y, preds[[nm]], quiet=TRUE)
    dl <- roc.test(rA, rX, method="delong")
    lr <- anova(fits$A, fits[[nm]], test="Chisq")
    out[[length(out)+1]] <- data.frame(
      outcome=lab, model=nm, n=nrow(dd), events=sum(y),
      auc_A=as.numeric(auc(rA)), auc_X=as.numeric(auc(rX)),
      dAUC=as.numeric(auc(rX))-as.numeric(auc(rA)),
      p_delong=dl$p.value,
      LR_chisq=lr$Deviance[2], LR_df=lr$Df[2], p_LR=lr$`Pr(>Chi)`[2],
      dAIC=AIC(fits[[nm]])-AIC(fits$A),
      NRI=nri_cont(preds$A, preds[[nm]], y),
      IDI=idi(preds$A, preds[[nm]], y), stringsAsFactors=FALSE)
  }
  # cluster bootstrap for dAUC of B3 and B4
  cl <- unique(dd$communityID); B <- 400
  store <- matrix(NA, B, 2, dimnames=list(NULL, c("B3","B4")))
  for (i in 1:B) {
    pick <- sample(cl, length(cl), replace=TRUE)
    idx <- unlist(lapply(pick, function(k) which(dd$communityID==k)))
    d2 <- dd[idx, ]; y2 <- y[idx]
    if (length(unique(y2))<2) next
    pA2 <- predict(glm(as.formula(paste("y ~",paste(c(COV,"g13_5"),collapse="+"))),data=d2,family=binomial()),type="response")
    for (nm in c("B3","B4")) {
      pX2 <- tryCatch(predict(glm(as.formula(paste("y ~",paste(c(COV,specs[[nm]]),collapse="+"))),data=d2,family=binomial()),type="response"), error=function(e) NULL)
      if (!is.null(pX2)) store[i,nm] <- as.numeric(auc(roc(y2,pX2,quiet=TRUE))-auc(roc(y2,pA2,quiet=TRUE)))
    }
  }
  cat(sprintf("\n=== %s (n=%d, events=%d) ===\n", lab, nrow(dd), sum(y)))
  for (nm in c("B3","B4")) {
    v <- store[,nm]; v <- v[is.finite(v)]
    cat(sprintf("  bootstrap dAUC %s: 95%%CI %+.4f to %+.4f (B=%d)\n", nm, quantile(v,.025), quantile(v,.975), length(v)))
  }
}
res <- do.call(rbind, out)
print(res, digits=4, row.names=FALSE)
write.csv(res, "out/S9_incremental_continuous.csv", row.names=FALSE)

sink()
