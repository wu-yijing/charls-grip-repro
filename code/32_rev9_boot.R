# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 32_rev9_boot.R -- cluster-bootstrap dAUC for every incremental model (B1-B5)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
m <- readRDS("out/jamda_full.rds")
sink("out/32_rev9_boot.txt", split = TRUE)
set.seed(20261024)

b <- get_base(m, "grp19")
b$g13_5 <- b$g13/5; b$g11_5 <- b$g11/5; b$chg_5 <- (b$g13-b$g11)/5
b$low11 <- b$low11_19; b$grp4 <- b$grp
d0 <- b[complete.cases(b[, c(COV,"g13_5","g11_5","chg_5","low11","grp4","oc")]), ]

specs <- list(B1=c("g13_5","low11"), B2=c("g13_5","grp4"),
              B3=c("g13_5","g11_5"), B4=c("g13_5","chg_5"),
              B5=c("g13_5","g11_5","low11"))
fA <- as.formula(paste("y ~", paste(c(COV,"g13_5"), collapse="+")))
B <- 300
for (outcome in c(1,2)) {
  dd <- d0[d0$oc %in% c(0,outcome), ]; y <- as.numeric(dd$oc==outcome)
  cl <- unique(dd$communityID)
  store <- matrix(NA, B, length(specs), dimnames=list(NULL, names(specs)))
  for (i in 1:B) {
    pick <- sample(cl, length(cl), replace=TRUE)
    idx <- unlist(lapply(pick, function(k) which(dd$communityID==k)))
    d2 <- dd[idx,]; y2 <- y[idx]
    if (length(unique(y2))<2) next
    pA <- predict(glm(fA, data=d2, family=binomial()), type="response")
    aA <- as.numeric(auc(roc(y2,pA,quiet=TRUE)))
    for (nm in names(specs)) {
      pX <- tryCatch(predict(glm(as.formula(paste("y ~",paste(c(COV,specs[[nm]]),collapse="+"))),
                                 data=d2, family=binomial()), type="response"), error=function(e) NULL)
      if (!is.null(pX)) store[i,nm] <- as.numeric(auc(roc(y2,pX,quiet=TRUE))) - aA
    }
  }
  cat("\n---", c("disability","death")[outcome], "---\n")
  for (nm in names(specs)) {
    v <- store[,nm]; v <- v[is.finite(v)]
    cat(sprintf("  %s: 95%%CI %+.4f to %+.4f (n_boot=%d)\n", nm, quantile(v,.025), quantile(v,.975), length(v)))
  }
  write.csv(store, sprintf("out/S9_boot_%s.csv", c("disability","death")[outcome]), row.names=FALSE)
}
cat("\ndone\n")

sink()
