source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages({library(nnet);library(dplyr)})
m <- readRDS("out/jamda_full.rds")
b0 <- get_base(m, "grp19")
f <- as.formula(paste("y ~", paste(c(COV,"grp"), collapse="+")))
run <- function(occol, lab){
  d <- b0; d$y <- factor(ifelse(d$dead18==1, 2, ifelse(!is.na(d[[occol]]), as.numeric(d[[occol]]==1), NA)),
                         levels=0:2, labels=c("alive_free","disabled","dead"))
  dc <- d[complete.cases(d[, c("y", COV, "grp")]), ]
  ff <- f; environment(ff) <- environment()   # f was built at top level; bind it to this frame
  fit <- multinom(ff, data=dc, trace=FALSE, maxit=500)
  cf <- summary(fit)$coefficients
  cat(sprintf("\n=== %s (n=%d) ===\n", lab, nrow(dc)))
  for (k in 1:2) for (g in c("remitted","incident","persistent"))
    cat(sprintf("%-9s %-11s %.2f\n", c("disabled","dead")[k], g, exp(cf[k, paste0("grp",g)])))
}
run("aD4", "disability requiring inability (ADL>=4)")
run("aC",  "ADL+IADL>=3 (11 items)")
# sex x group interaction
dc <- b0[complete.cases(b0[, c("y", COV, "grp")]), ]
gd <- model.matrix(~ grp - 1, dc)
colnames(gd) <- sub("grp","",colnames(gd))
gd <- gd[, c("remitted","incident","persistent")]
X0 <- cbind(1, as.matrix(dc[,COV]), gd)
X1 <- cbind(X0, sweep(gd, 1, dc$male, "*"))
colnames(X1)[(ncol(X0)+1):ncol(X1)] <- paste0(colnames(gd), "_x_male")
m0 <- multinom(y ~ ., data=data.frame(y=dc$y, X0), trace=FALSE, maxit=500)
m1 <- multinom(y ~ ., data=data.frame(y=dc$y, X1), trace=FALSE, maxit=500)
ll0 <- as.numeric(logLik(m0)); ll1 <- as.numeric(logLik(m1))
lr <- 2*(ll1-ll0); df <- ncol(X1)-ncol(X0)
cat(sprintf("\nsex x group interaction: LR chi2(%d) = %.2f, p = %.3f\n", df, lr, pchisq(lr, df, lower.tail=FALSE)))
