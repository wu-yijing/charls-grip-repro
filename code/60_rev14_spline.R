# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 60_rev14_spline.R ------------------------------------------------------------------
# Revision-14 analyses requested by the reviewer report:
#   (A) reproduction of the primary (row 1) and row 17 estimates
#   (B) current value entered as a natural cubic spline (3 and 4 df) and as quintiles
#   (C) ordinal term (number of low readings) with pooled CI and test, for a new S5 row
#   (D) Brenner matrix correction for misclassification using the latent-class Se/Sp
# Output: out/60_rev14_spline.txt
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(mice); library(nnet); library(dplyr); library(splines)})
m <- readRDS("out/jamda_full.rds")
sink("out/60_rev14_spline.txt", split = TRUE)

mn_vcov <- function(beta, X, y, w, cluster){
  K <- 2; p <- ncol(X); N <- nrow(X)
  eta <- X %*% t(beta); P <- exp(cbind(0, eta)); P <- P/rowSums(P)
  Y <- sapply(1:2, function(k) as.numeric(y==k))
  A <- matrix(0, K*p, K*p)
  Wm <- list(P[,2]*(1-P[,2]), P[,3]*(1-P[,3]), -P[,2]*P[,3])
  bl <- matrix(c(1,3,3,2), 2, 2, byrow=TRUE)
  for (k in 1:K) for (l in 1:K) {
    ww <- w * Wm[[bl[k,l]]]
    A[((k-1)*p+1):(k*p), ((l-1)*p+1):(l*p)] <- crossprod(X, X*ww)
  }
  S <- matrix(0, N, K*p)
  for (k in 1:K) S[, ((k-1)*p+1):(k*p)] <- (w*(Y[,k]-P[,k+1])) * X
  Sc <- rowsum(S, group=cluster, reorder=FALSE); G <- nrow(Sc)
  Ai <- tryCatch(solve(A), error=function(e) MASS::ginv(A))
  list(V = (G/(G-1))*((N-1)/(N-K*p)) * Ai %*% crossprod(Sc) %*% Ai)
}

rubin <- function(bs, vs){
  mm <- length(bs); b <- Reduce(`+`, bs)/mm
  W <- Reduce(`+`, vs)/mm
  B <- Reduce(`+`, lapply(bs, function(z) (z-b)^2))/(mm-1)
  T <- W + (1+1/mm)*B
  lam <- (1+1/mm)*B/pmax(T,1e-12)
  df <- (mm-1)/pmax(lam,1e-12)^2
  list(b=b, se=sqrt(T), df=ifelse(is.finite(df) & df>0, df, 1e5))
}

b0 <- get_base(m, "grp19")
d0 <- b0 %>% filter(!is.na(oc), !is.na(w13), w13 > 0) %>%
  mutate(nr = low11_19 + low13_19, nr_low = nr, edu_f = factor(edu),
         wsc = w13/mean(w13, na.rm=TRUE), g13_5 = g13/5)
d0$nr_low <- d0$nr_low - 0

cat("n =", nrow(d0), " clusters =", length(unique(d0$communityID)), "\n")
cat("deaths =", sum(d0$oc==2), " disabled =", sum(d0$oc==1), "\n\n")

## ---- fixed spline knots (deterministic across imputations) ----
q <- quantile(d0$g13, c(.05,.10,.35,.50,.65,.90,.95))
bc <- c(min(d0$g13), max(d0$g13))
cat("g13 quantiles: "); print(round(q,1)); cat("boundary:", bc, "\n\n")
B3 <- ns(d0$g13, knots=as.numeric(q[c("35%","65%")]), Boundary.knots=bc)   # 3 df
B4 <- ns(d0$g13, knots=as.numeric(q[c("10%","50%","90%")]), Boundary.knots=bc) # 4 df
colnames(B3) <- paste0("spl", 1:ncol(B3))
colnames(B4) <- paste0("spl", 1:ncol(B4))

imp <- mice(d0[, c("age10","male","bmi","nchronic","cesd","smoke","edu_f","grp","y","nr_low")],
            m=20, seed=20261023, printFlag=FALSE)

pool <- function(extra, type){
  bs <- list(); vs <- list(); cols <- NULL
  for (i in 1:20) {
    di <- complete(imp, i)
    di$y   <- factor(di$y, levels=c("alive_free","disabled","dead"))
    di$grp <- relevel(factor(di$grp, levels=GRP), ref="never")
    di$edu <- as.numeric(as.character(di$edu_f))
    for (kk in 1:4) di[[paste0("edu",kk)]] <- as.numeric(di$edu==kk)
    di$communityID <- d0$communityID; di$wsc <- d0$wsc
    di$g13 <- d0$g13; di$g13_5 <- d0$g13/5
    if (type=="spl3") di <- cbind(di, as.data.frame(B3)); 
    if (type=="spl4") di <- cbind(di, as.data.frame(B4))
    X <- cbind(1, as.matrix(di[, COV]))
    if (type != "ord") for (g in c("remitted","incident","persistent")) X <- cbind(X, as.numeric(di$grp==g))
    add <- switch(type,
      none  = character(0),
      lin   = "g13_5",
      spl3  = colnames(as.data.frame(B3)),
      spl4  = colnames(as.data.frame(B4)),
      quin  = paste0("q",2:5),
      ord   = "nr_low")
    if (type=="quin") for (kk in 1:4) di[[paste0("q",kk+1)]] <- as.numeric(as.integer(cut(di$g13, breaks=quantile(d0$g13, seq(0,1,.2)), include.lowest=TRUE, labels=FALSE))==kk+1)
    if (length(add)) X <- cbind(X, as.matrix(di[, add]))
    cn <- if (type=="ord") c("const", COV, add) else c("const", COV, paste0("grp", c("remitted","incident","persistent")), add)
    colnames(X) <- cn
    rhs <- if (type=="ord") c(COV, "nr_low") else c(COV, "grp", add)
    f <- as.formula(paste("y ~", paste(rhs, collapse="+")))
    fit <- multinom(f, data=di, weights=di$wsc, trace=FALSE, maxit=500)
    bb <- matrix(coef(fit), nrow=2, byrow=FALSE); colnames(bb) <- cn
    vc <- mn_vcov(bb, X, as.integer(di$y)-1, di$wsc, di$communityID)
    bs[[i]] <- as.numeric(t(bb)); vs[[i]] <- as.numeric(diag(vc$V))
    cols <- cn
  }
  pr <- rubin(bs, vs); p <- length(cols)
  list(B=matrix(pr$b,2,byrow=TRUE), S=matrix(pr$se,2,byrow=TRUE),
       D=matrix(pr$df,2,byrow=TRUE), cols=cols, p=p)
}
getc <- function(r, nm){ j <- which(r$cols==nm)[1]
  b <- r$B[,j]; s <- r$S[,j]; dd <- r$D[,j]
  data.frame(outcome=c("disabled","dead"), RRR=exp(b),
             lo=exp(b-qt(.975,dd)*s), hi=exp(b+qt(.975,dd)*s),
             p=2*pt(-abs(b/s),dd), row.names=NULL) }

cat("#################### A. REPRODUCTION ####################\n")
r0 <- pool(NULL, "none")
cat("primary (covariates, no current value):\n"); print(getc(r0,"grppersistent"), digits=3)
r17 <- pool(NULL, "lin")
cat("\nrow 17 (covariates + 2013 grip per 5 kg, linear):\n"); print(getc(r17,"grppersistent"), digits=3)

cat("\n#################### B. CURRENT VALUE AS SPLINE ####################\n")
r3 <- pool(NULL, "spl3"); cat("3 df natural spline:\n"); print(getc(r3,"grppersistent"), digits=3)
r4 <- pool(NULL, "spl4"); cat("\n4 df natural spline:\n"); print(getc(r4,"grppersistent"), digits=3)
rq <- pool(NULL, "quin"); cat("\nquintiles:\n"); print(getc(rq,"grppersistent"), digits=3)

cat("\n#################### C. ORDINAL TERM ####################\n")
ro <- pool(NULL, "ord")
cat("per additional low reading (number of low readings, 0-2):\n"); print(getc(ro,"nr_low"), digits=3)

## unweighted, complete-case ordinal fit for the likelihood-ratio test (df check)
dc <- d0[complete.cases(d0[,c("y",COV,"grp","nr_low")]),]
COVl <- COV
m4 <- multinom(as.formula(paste("y ~", paste(c(COVl,"grp"), collapse="+"))), data=dc, weights=dc$wsc, trace=FALSE, maxit=500)
m2 <- multinom(as.formula(paste("y ~", paste(c(COVl,"nr_low"), collapse="+"))), data=dc, weights=dc$wsc, trace=FALSE, maxit=500)
lr <- 2*(as.numeric(logLik(m4)) - as.numeric(logLik(m2)))
cat(sprintf("\nLR chi2 = %.3f | df from parameter count = %d | P = %.4f\n", lr,
            attr(logLik(m4),"df") - attr(logLik(m2),"df"), pchisq(lr, attr(logLik(m4),"df") - attr(logLik(m2),"df"), lower.tail=FALSE)))
cat(sprintf("(df if a single-outcome contrast was intended: %d, P = %.4f)\n",
            (attr(logLik(m4),"df") - attr(logLik(m2),"df"))/2, pchisq(lr, (attr(logLik(m4),"df")-attr(logLik(m2),"df"))/2, lower.tail=FALSE)))
cat("n in LR test:", nrow(dc), "\n")

cat("\n#################### D. BRENNER MATRIX CORRECTION ####################\n")
Se <- 0.594; Sp <- 0.950; SeCI <- c(0.540, 0.637); SpCI <- c(0.943, 0.961)
corr <- function(n_ev_low, n_ev_tot, n_ne_low, n_ne_tot, se, sp){
  f <- function(n1, N){ (n1 - N*(1-sp))/(se+sp-1) }
  m1 <- f(n_ev_low, n_ev_tot); m0 <- f(n_ne_low, n_ne_tot)
  c(rr=(m1/n_ev_tot)/(m0/n_ne_tot), miss1=(m1), miss0=(m0))
}
tab <- function(expo, ev){
  data.frame(
    ev_low = sum(expo & ev),  ev_tot = sum(ev),
    ne_low = sum(expo & !ev), ne_tot = sum(!ev))
}
out <- list()
for (oc_name in c("dead","disabled")) {
  ev <- d0$oc == ifelse(oc_name=="dead", 2, 1)
  # exposure 1: low in 2013 (single reading)
  t1 <- tab(d0$low13_19==1, ev)
  r_a <- corr(t1$ev_low, t1$ev_tot, t1$ne_low, t1$ne_tot, Se, Sp)
  r_lo <- corr(t1$ev_low, t1$ev_tot, t1$ne_low, t1$ne_tot, SeCI[1], SpCI[2])
  r_hi <- corr(t1$ev_low, t1$ev_tot, t1$ne_low, t1$ne_tot, SeCI[2], SpCI[1])
  # exposure 2: persistent vs never-low (Se^2, Sp^2)
  sub <- d0$grp19 %in% c("never","persistent")
  t2 <- tab(d0$grp19=="persistent" & sub, ev & sub)
  r2  <- corr(t2$ev_low, t2$ev_tot, t2$ne_low, t2$ne_tot, Se^2, Sp^2)
  r2lo<- corr(t2$ev_low, t2$ev_tot, t2$ne_low, t2$ne_tot, SeCI[1]^2, SpCI[2]^2)
  r2hi<- corr(t2$ev_low, t2$ev_tot, t2$ne_low, t2$ne_tot, SeCI[2]^2, SpCI[1]^2)
  # uncorrected
  u1 <- (t1$ev_low/t1$ev_tot)/(t1$ne_low/t1$ne_tot)
  u2 <- (t2$ev_low/t2$ev_tot)/(t2$ne_low/t2$ne_tot)
  out[[oc_name]] <- data.frame(outcome=oc_name,
    single_obs=u1, single_cor=r_a["rr"], single_lo=r_lo["rr"], single_hi=r_hi["rr"],
    pers_obs=u2, pers_cor=r2["rr"], pers_lo=r2lo["rr"], pers_hi=r2hi["rr"])
  cat("\n--- ", oc_name, " ---\n", sep="")
  cat(sprintf("single reading (2013 low vs not): observed RR %.2f | corrected %.2f (range %.2f-%.2f)\n", u1, r_a["rr"], r_lo["rr"], r_hi["rr"]))
  cat(sprintf("persistent vs never (Se^2=%.3f, Sp^2=%.3f): observed RR %.2f | corrected %.2f (range %.2f-%.2f)\n",
              Se^2, Sp^2, u2, r2["rr"], r2lo["rr"], r2hi["rr"]))
  cat(sprintf("  implied true low counts: single %.0f of %.0f exposed (observed %d)\n", r_a["miss1"], t1$ev_tot, t1$ev_low))
}
print(do.call(rbind, out), digits=2, row.names=FALSE)
cat("\nNote: corrected counts can be negative when the observed prevalence is below 1-Sp; such cells are flagged above.\n")

sink()
