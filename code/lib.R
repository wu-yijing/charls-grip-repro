# lib.R -- shared helpers for the JAMDA re-analyses
suppressMessages({library(nnet); library(sandwich); library(dplyr); library(pROC)})

COV <- c("age10","male","bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4")
GRP <- c("never","remitted","incident","persistent")
CN  <- c(never="\u59cb\u7ec8\u6b63\u5e38", remitted="\u5df2\u56de\u5347",
         incident="\u65b0\u53d1\u4f4e\u63e1\u529b", persistent="\u6301\u7eed\u4f4e\u63e1\u529b")

get_base <- function(m, grpcol){
  b <- m %>% filter(a13 %in% 0, !is.na(.data[[grpcol]]), !is.na(age13), !is.na(male))
  b$grp <- relevel(factor(b[[grpcol]], levels=GRP), ref="never")
  b$y   <- factor(b$oc, levels=0:2, labels=c("alive_free","disabled","dead"))
  b
}

# multinomial logit; returns tidy RRR table
fit_mn <- function(d, covars, cluster=NULL){
  d <- d[complete.cases(d[, c("y", covars, "grp")]), ]
  f <- as.formula(paste("y ~", paste(c(covars,"grp"), collapse="+")))
  fit <- multinom(f, data=d, trace=FALSE, maxit=200)
  cf <- summary(fit)$coefficients; se <- summary(fit)$standard.errors
  if (is.null(dim(cf))) { cf <- matrix(cf, nrow=1); se <- matrix(se, nrow=1) }
  V <- NULL
  if (!is.null(cluster)) {
    V <- tryCatch(vcovCL(fit, cluster=cluster[complete.cases(d[, c("y", covars, "grp")])], type="HC0"),
                  error=function(e) NULL)
  }
  out <- list(fit=fit, coef=cf, se=se, V=V, n=nrow(d), data=d)
  out
}

rrr_rows <- function(res, label){
  cf <- res$coef; se <- res$se
  eqs <- rownames(cf)
  rows <- list()
  for (eq in seq_len(nrow(cf))) {
    em <- c("alive_free","disabled","dead")[match(eqs[eq], c("alive_free","disabled","dead"))]
    lab <- if (eqs[eq]=="disabled") "disabled" else if (eqs[eq]=="dead") "dead" else eqs[eq]
    for (g in c("remitted","incident","persistent")) {
      k <- paste0("grp", g)
      if (!(k %in% colnames(cf))) next
      b <- cf[eq,k]; s <- se[eq,k]
      rows[[length(rows)+1]] <- data.frame(model=label, outcome=lab, group=g,
        RRR=exp(b), lo=exp(b-1.96*s), hi=exp(b+1.96*s),
        p=2*pnorm(-abs(b/s)), n=res$n, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}

fmt_rrr <- function(r, digits=2)
  sprintf("%.2f (%.2f-%.2f)", r$RRR, r$lo, r$hi)

# predicted absolute risks at given ages, reference profile = mean covariates, 50% male
pred_risk <- function(res, ages=c(60,65,70,75,80), male=0.5){
  d <- res$data; cf <- res$coef
  cols <- colnames(cf)
  out <- list()
  for (a in ages) {
    for (g in GRP) {
      x <- setNames(rep(0, length(cols)), cols)
      x["(Intercept)"] <- 1
      x["male"] <- male
      x["age10"] <- (a-65)/10
      for (cc in c("bmi","nchronic","cesd","smoke","edu1","edu2","edu3","edu4"))
        if (cc %in% cols) x[cc] <- mean(d[[cc]], na.rm=TRUE)
      if (g != "never" && paste0("grp",g) %in% cols) x[paste0("grp",g)] <- 1
      eta <- as.numeric(cf %*% x)                 # length 2: disabled, dead
      p <- exp(c(0, eta)); p <- p/sum(p)
      out[[length(out)+1]] <- data.frame(age=a, group=g,
        disabled=100*p[2], dead=100*p[3], stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, out)
}

# 2x2 transition table + kappa
trans_table <- function(m, l11, l13){
  d <- m[!is.na(m[[l11]]) & !is.na(m[[l13]]), ]
  tab <- table(factor(d[[l11]], 0:1), factor(d[[l13]], 0:1))
  po <- sum(diag(tab))/sum(tab)
  pe <- sum(rowSums(tab)*colSums(tab))/sum(tab)^2
  list(tab=tab, po=po, pe=pe, kappa=(po-pe)/(1-pe), n=nrow(d))
}
