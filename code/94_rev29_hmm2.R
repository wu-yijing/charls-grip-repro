# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
## 94_rev29_hmm2.R -- two-state latent Markov model fitted directly to the 2x2x2 table
suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male)) %>%
  mutate(cut=ifelse(male==1,28,18))
b$l11 <- as.numeric(b$g11<b$cut); b$l13 <- as.numeric(b$g13<b$cut); b$l15 <- as.numeric(b$g15<b$cut)
cc <- b %>% filter(!is.na(l11), !is.na(l13), !is.na(l15))
cat("complete three-wave n =", nrow(cc), "\n")
cell <- cc %>% group_by(l11,l13,l15) %>% summarise(n=n(), .groups="drop")
print(as.data.frame(cell))
N <- sum(cell$n)
## index: y1*4 + y2*2 + y3 + 1  (0=normal,1=low)
cnt <- numeric(8)
for (i in seq_len(nrow(cell))) cnt[cell$l11[i]*4 + cell$l13[i]*2 + cell$l15[i] + 1] <- cell$n[i]

prob <- function(par, tv = FALSE) {
  pi1 <- plogis(par[1])
  if (!tv) { a01 <- plogis(par[2]); a10 <- plogis(par[3])
             A1 <- matrix(c(1-a01, a01, a10, 1-a10), 2, byrow=TRUE); A2 <- A1
             e1 <- plogis(par[4]); e0 <- plogis(par[5]); k <- 5
  } else {   a01a <- plogis(par[2]); a10a <- plogis(par[3]); a01b <- plogis(par[4]); a10b <- plogis(par[5])
             A1 <- matrix(c(1-a01a, a01a, a10a, 1-a10a), 2, byrow=TRUE)
             A2 <- matrix(c(1-a01b, a01b, a10b, 1-a10b), 2, byrow=TRUE)
             e1 <- plogis(par[6]); e0 <- plogis(par[7]); k <- 7 }
  E <- cbind(c(1-e0, 1-e1), c(e0, e1))            # rows: states (1 normal, 2 low), cols: y=0/1
  pin <- c(1-pi1, pi1)
  p <- numeric(8)
  for (y1 in 0:1) for (y2 in 0:1) for (y3 in 0:1) {
    s <- 0
    for (s1 in 1:2) for (s2 in 1:2) for (s3 in 1:2)
      s <- s + pin[s1]*A1[s1,s2]*A2[s2,s3]*
                E[s1, y1+1]*E[s2, y2+1]*E[s3, y3+1]
    p[y1*4+y2*2+y3+1] <- s
  }
  list(p=p, A1=A1, A2=A2, e1=e1, e0=e0, pi1=pi1, k=k)
}
nll <- function(par, tv) { p <- prob(par, tv)$p; if (any(p<=0)) return(1e9); -sum(cnt*log(p)) }

best <- NULL
for (s in 1:60) {
  set.seed(500+s)
  st <- c(qlogis(runif(1,.03,.25)), qlogis(runif(1,.0005,.05)), qlogis(runif(1,.02,.3)),
          qlogis(runif(1,.5,.95)), qlogis(runif(1,.005,.12)))
  o <- try(optim(st, nll, tv=FALSE, method="BFGS", control=list(maxit=2000, reltol=1e-12)), silent=TRUE)
  if (inherits(o,"try-error")) next
  if (is.null(best) || o$value < best$value) best <- o
}
cat(sprintf("\nstationary transitions: logL = %.4f (k=5, df=%d), converged %d\n", -best$value, 7-5, best$convergence))
P <- prob(best$par, FALSE)
cat(sprintf("pi(low) %.4f ; P(normal->low) %.5f ; P(low->normal) %.4f ; sens %.4f ; 1-spec %.4f\n",
            P$pi1, P$A1[1,2], P$A1[2,1], P$e1, P$e0))
cat("fitted cell proportions vs observed:/n")
print(round(data.frame(fitted=as.numeric(P$p), observed=as.numeric(cnt)/N), 4))
lr <- 2*(sum(ifelse(cnt>0, cnt*log(cnt/N), 0)) - (-best$value))
cat(sprintf("GOF LR = %.2f on %d df, p = %.4f\n", lr, 2, pchisq(lr, 2, lower.tail=FALSE)))

## time-varying transitions (just identified, k = 7)
best2 <- NULL
for (s in 1:60) {
  set.seed(900+s)
  st <- c(qlogis(runif(1,.03,.25)), qlogis(runif(1,.0005,.05)), qlogis(runif(1,.02,.3)),
          qlogis(runif(1,.0005,.08)), qlogis(runif(1,.02,.4)),
          qlogis(runif(1,.5,.95)), qlogis(runif(1,.005,.12)))
  o <- try(optim(st, nll, tv=TRUE, method="BFGS", control=list(maxlist=2000, maxit=4000, reltol=1e-12)), silent=TRUE)
  if (inherits(o,"try-error")) o <- try(optim(st, nll, tv=TRUE, method="Nelder-Mead", control=list(maxit=20000, reltol=1e-12)), silent=TRUE)
  if (inherits(o,"try-error")) next
  if (is.null(best2) || o$value < best2$value) best2 <- o
}
if (!is.null(best2)) {
  Q <- prob(best2$par, TRUE)
  cat(sprintf("\ntime-varying transitions: logL = %.4f (k=7, df=0)\n", -best2$value))
  cat(sprintf("wave1->2: P(normal->low) %.5f ; P(low->normal) %.4f\n", Q$A1[1,2], Q$A1[2,1]))
  cat(sprintf("wave2->3: P(normal->low) %.5f ; P(low->normal) %.4f\n", Q$A2[1,2], Q$A2[2,1]))
  cat(sprintf("sens %.4f ; 1-spec %.4f ; pi(low) %.4f\n", Q$e1, Q$e0, Q$pi1))
}
obs_rev <- sum(cnt[5:8]) ; rev11_13 <- (cnt[5]+cnt[7])/obs_rev
cat(sprintf("\nobserved reversion 11->13 (complete cases) = %.4f\n", rev11_13))
cat(sprintf("observed reversion 13->15 = %.4f\n", (cnt[5]+cnt[6])/sum(cnt[3:6])))
