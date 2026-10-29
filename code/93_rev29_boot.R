# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
## 93_rev29_boot.R -- cluster bootstrap (communities) for the agreement indices in Table 1
suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male), !is.na(communityID)) %>%
  mutate(cut=ifelse(male==1,28,18), l11=as.numeric(g11<cut), l13=as.numeric(g13<cut))
cat("n =", nrow(b), " communities =", length(unique(b$communityID)), "\n")
a <- sum(b$l11==1 & b$l13==1); bb <- sum(b$l11==1 & b$l13==0)
cc <- sum(b$l11==0 & b$l13==1); d <- sum(b$l11==0 & b$l13==0)
idx <- function(a,bb,cc,d){
  n <- a+bb+cc+d; po <- (a+d)/n
  pe <- ((a+bb)*(a+cc) + (cc+d)*(bb+d))/n^2
  kap <- (po-pe)/(1-pe)
  c(agree=po, kappa=kap, pabak=2*po-1, pa=2*a/(2*a+bb+cc), na=2*d/(2*d+bb+cc),
    prev=abs(a-d)/n, bias=abs(bb-cc)/n, rev=bb/(a+bb), inc=cc/(cc+d), a=a, b=bb, c=cc, d=d)
}
pt <- idx(a,bb,cc,d)
cat("point:", sprintf("%s=%.4f", names(pt), pt), "\n\n")
set.seed(20261028)
com <- split(seq_len(nrow(b)), b$communityID); ncom <- length(com)
R <- 1000; out <- matrix(NA, R, length(pt)); colnames(out) <- names(pt)
for (r in 1:R) {
  pick <- sample(ncom, ncom, replace=TRUE)
  ii <- unlist(com[pick], use.names=FALSE)
  bb2 <- b[ii, ]
  out[r, ] <- idx(sum(bb2$l11==1 & bb2$l13==1), sum(bb2$l11==1 & bb2$l13==0),
                  sum(bb2$l11==0 & bb2$l13==1), sum(bb2$l11==0 & bb2$l13==0))
}
ci <- apply(out, 2, quantile, c(.025,.975), na.rm=TRUE)
for (k in c("agree","kappa","pabak","pa","na","rev","inc")) {
  cat(sprintf("%-6s %.4f  cluster bootstrap 95%% CI %.4f-%.4f\n", k, pt[k], ci[1,k], ci[2,k]))
}
## for comparison: individual-level bootstrap
out2 <- matrix(NA, 500, length(pt)); colnames(out2) <- names(pt)
for (r in 1:500) {
  ii <- sample(nrow(b), nrow(b), replace=TRUE); bb2 <- b[ii, ]
  out2[r,] <- idx(sum(bb2$l11==1 & bb2$l13==1), sum(bb2$l11==1 & bb2$l13==0),
                  sum(bb2$l11==0 & bb2$l13==1), sum(bb2$l11==0 & bb2$l13==0))
}
ci2 <- apply(out2, 2, quantile, c(.025,.975), na.rm=TRUE)
cat("\nindividual-level bootstrap (for the design-effect check):\n")
for (k in c("kappa","pa","na","rev")) cat(sprintf("%-6s %.4f-%.4f\n", k, ci2[1,k], ci2[2,k]))
cat(sprintf("\nwidth ratio cluster/individual: kappa %.2f  pa %.2f  na %.2f  rev %.2f\n",
            (ci[2,"kappa"]-ci[1,"kappa"])/(ci2[2,"kappa"]-ci2[1,"kappa"]),
            (ci[2,"pa"]-ci[1,"pa"])/(ci2[2,"pa"]-ci2[1,"pa"]),
            (ci[2,"na"]-ci[1,"na"])/(ci2[2,"na"]-ci2[1,"na"]),
            (ci[2,"rev"]-ci[1,"rev"])/(ci2[2,"rev"]-ci2[1,"rev"])))
