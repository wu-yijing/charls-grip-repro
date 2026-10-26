source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages({library(mice);library(nnet);library(dplyr)})
m <- readRDS("out/jamda_full.rds")
b0 <- get_base(m, "grp19")
d0 <- b0 %>% filter(!is.na(oc)) %>% mutate(edu_f = factor(edu))
# (a) 完整病例、未加权、聚类稳健
dc <- d0[complete.cases(d0[, c("y", COV, "grp")]), ]
X <- cbind(1, as.matrix(dc[, COV]), as.numeric(dc$grp=="remitted"),
           as.numeric(dc$grp=="incident"), as.numeric(dc$grp=="persistent"))
colnames(X) <- c("const", COV, paste0("grp", c("remitted","incident","persistent")))
f <- as.formula(paste("y ~", paste(c(COV,"grp"), collapse="+")))
fit <- multinom(f, data=dc, trace=FALSE, maxit=500)
b <- matrix(coef(fit), 2, byrow=FALSE); colnames(b) <- colnames(X)
cat("=== (a) 未加权完整病例 n=", nrow(dc), " ===\n")
for (k in 1:2) for (j in 12:14)
  cat(sprintf("%-9s %-11s %.2f\n", c("disabled","dead")[k], c("remitted","incident","persistent")[j-11], exp(b[k,j])))
# (b) 未加权 MI
imp <- mice(d0[, c(COV,"edu_f","grp","y")], m=20, seed=20261023, printFlag=FALSE)
bs <- list(); vs <- list()
for (i in 1:20) {
  di <- complete(imp, i); di$y <- factor(di$y, levels=c("alive_free","disabled","dead"))
  di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
  di$edu <- as.numeric(as.character(di$edu_f)); for (k in 1:4) di[[paste0("edu",k)]] <- as.numeric(di$edu==k)
  fit <- multinom(f, data=di, weights=rep(1, nrow(di)), trace=FALSE, maxit=500)
  cf <- summary(fit)$coefficients; se <- summary(fit)$standard.errors
  bs[[i]] <- as.numeric(t(cf)); vs[[i]] <- as.numeric(t(se^2))
}
b2 <- Reduce(`+`, bs)/20; W <- Reduce(`+`, vs)/20
Bv <- Reduce(`+`, lapply(bs, function(z) (z-b2)^2))/19; T <- W + (1+1/20)*Bv
p <- length(b2)/2
cat("\n=== (b) 未加权 MI（20 套）===\n")
for (k in 1:2) for (j in 12:14) {
  idx <- (k-1)*p + j
  cat(sprintf("%-9s %-11s %.2f (%.2f-%.2f)\n", c("disabled","dead")[k],
      c("remitted","incident","persistent")[j-11], exp(b2[idx]),
      exp(b2[idx]-1.96*sqrt(T[idx])), exp(b2[idx]+1.96*sqrt(T[idx]))))
}
