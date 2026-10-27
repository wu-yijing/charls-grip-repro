# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R) and, where noted, the synthetic
# raw files under data/raw_synthetic/. Nothing here reads CHARLS data.
# 75_kappa_cluster.R -- kappa CI: individual versus community bootstrap
#
# The archived script 21_enhanced_A.R draws `sample.int(n, n, TRUE)`, i.e. it resamples
# individuals. Every other interval in the paper is clustered on community, because grip
# and the outcome are correlated within a community. This script reports both, so the
# manuscript can state which unit the kappa interval uses, and switch if the difference matters.
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(dplyr)})
m <- readRDS("out/jamda_full.rds")
sink("out/75_kappa_cluster.txt", split = TRUE)

b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
cat("analytic sample n =", nrow(b), "| communities =", length(unique(b$communityID)), "\n")

l1 <- as.integer(b$low11_19); l3 <- as.integer(b$low13_19)
ok <- !is.na(l1) & !is.na(l3)
l1 <- l1[ok]; l3 <- l3[ok]; cid <- b$communityID[ok]
n <- length(l1); cat("classifiable pairs =", n, "\n")

kap <- function(a, bb, cc, dd){
  N <- a + bb + cc + dd
  po <- (a + dd)/N
  pe <- ((a + bb)*(a + cc) + (cc + dd)*(bb + dd))/N^2
  (po - pe)/(1 - pe)
}
tab <- table(factor(l1, 0:1), factor(l3, 0:1))
a <- tab[2,2]; bb <- tab[2,1]; cc <- tab[1,2]; dd <- tab[1,1]
k0 <- kap(a, bb, cc, dd)
cat(sprintf("\n2x2: a=%d b=%d c=%d d=%d\n", a, bb, cc, dd))
cat(sprintf("kappa = %.3f | observed agreement %.1f%%\n", k0, 100*(a+dd)/n))

set.seed(20261027)
B <- 2000

## ---- (1) individual bootstrap, as in 21_enhanced_A.R ----
cell <- l1*2 + l3 + 1L
ks <- numeric(B)
for (i in 1:B) {
  s  <- sample.int(n, n, TRUE)
  t2 <- matrix(tabulate(cell[s], nbins = 4), 2, 2)
  ks[i] <- kap(t2[2,2], t2[2,1], t2[1,2], t2[1,1])
}
ci_i <- quantile(ks, c(.025, .975), na.rm = TRUE)
cat(sprintf("\nindividual bootstrap : %.3f (%.3f-%.3f)  SE %.4f\n", k0, ci_i[1], ci_i[2], sd(ks)))

## ---- (2) community cluster bootstrap ----
byc <- split(seq_len(n), cid)
ids <- names(byc)
kc <- numeric(B)
for (i in 1:B) {
  pick <- sample(ids, length(ids), replace = TRUE)
  idx  <- unlist(byc[pick], use.names = FALSE)
  t2 <- table(factor(l1[idx], 0:1), factor(l3[idx], 0:1))
  kc[i] <- kap(t2[2,2], t2[2,1], t2[1,2], t2[1,1])
}
ci_c <- quantile(kc, c(.025, .975), na.rm = TRUE)
cat(sprintf("community bootstrap: %.3f (%.3f-%.3f)  SE %.4f\n", k0, ci_c[1], ci_c[2], sd(kc)))
cat(sprintf("\nvariance ratio cluster/individual = %.2f\n", (sd(kc)/sd(ks))^2))

## ---- (3) same two units for the second interval (2013-2015) ----
b2 <- b[!is.na(b$g15), ]
b2$low15_19 <- as.integer(b2$g15 < ifelse(b2$male == 1, 28, 18))
l1b <- as.integer(b2$low13_19); l3b <- as.integer(b2$low15_19); cidb <- b2$communityID
okb <- !is.na(l1b) & !is.na(l3b); l1b <- l1b[okb]; l3b <- l3b[okb]; cidb <- cidb[okb]
nb <- length(l1b)
set.seed(20261028)
bycb <- split(seq_len(nb), cidb); idsb <- names(bycb)
kb <- numeric(1000)
tab2 <- table(factor(l1b, 0:1), factor(l3b, 0:1))
k0b <- kap(tab2[2,2], tab2[2,1], tab2[1,2], tab2[1,1])
for (i in 1:1000) {
  pick <- sample(idsb, length(idsb), replace = TRUE)
  idx  <- unlist(bycb[pick], use.names = FALSE)
  t2 <- table(factor(l1b[idx], 0:1), factor(l3b[idx], 0:1))
  kb[i] <- kap(t2[2,2], t2[2,1], t2[1,2], t2[1,1])
}
cat(sprintf("\n2013-2015: n = %d | kappa %.3f, community bootstrap (%.3f-%.3f)\n",
            nb, k0b, quantile(kb,.025)[[1]], quantile(kb,.975)[[1]]))
