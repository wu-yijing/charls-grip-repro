# Ported for the reproducibility package: run from the repository root.
# Inputs: out/jamda_full.rds (built by 01_build.R). Nothing here reads CHARLS data.
# 95_rev31_age.R —— AWGS 2019 切点的年龄适配性：分层一致性指标（社区整群 bootstrap）
suppressMessages(library(dplyr))
set.seed(20261029)

m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(a13 == 0, !is.na(grp19), !is.na(age13), !is.na(male)) %>%
  mutate(band = ifelse(age13 < 50, "45-49", ifelse(age13 < 65, "50-64", "65+")))

calc <- function(l11, l13, dat) {
  t  <- table(factor(dat[[l11]], 0:1), factor(dat[[l13]], 0:1))
  LL <- t[2,2]; LN <- t[2,1]; NL <- t[1,2]; NN <- t[1,1]
  n  <- sum(t); po <- (LL + NN)/n
  pe <- (sum(t[2,]) * sum(t[,2]) + sum(t[1,]) * sum(t[,1]))/n^2
  c(n=n, LL=LL, LN=LN, NL=NL, NN=NN, po=po, kappa=(po-pe)/(1-pe),
    pa=2*LL/(2*LL+LN+NL), na=2*NN/(2*NN+LN+NL), pabak=2*po-1,
    prev11=sum(t[2,])/n, rev=LN/(LL+LN), inc=NL/(NL+NN))
}

boot <- function(l11, l13, dat, B=400) {
  parts <- split(seq_len(nrow(dat)), dat$communityID)
  keys  <- names(parts); out <- matrix(NA_real_, B, 3)
  for (i in seq_len(B)) {
    rows <- unlist(parts[sample(keys, length(keys), TRUE)], use.names=FALSE)
    v <- calc(l11, l13, dat[rows, ])
    out[i, ] <- c(v["kappa"], v["rev"], v["pa"])
  }
  apply(out, 2, quantile, c(.025,.975), na.rm=TRUE)
}

lab <- function(x) sprintf("%.1f", 100*x)
ci  <- function(q) sprintf("%.1f-%.1f", 100*q[1], 100*q[2])

cat("n =", nrow(b), " | <65:", lab(mean(b$age13<65)), "% | 45-49:", sum(b$band=="45-49"),
    " 50-64:", sum(b$band=="50-64"), " 65+:", sum(b$band=="65+"),
    " | age13<45:", sum(b$age13<45), " age11<45:", sum(b$age11<45, na.rm=TRUE), "\n\n")

specs <- list(
  list(tag="all (2019)",   dat=b,                    l11="low11_19", l13="low13_19"),
  list(tag="age <65 (2019)", dat=b[b$age13<65, ],    l11="low11_19", l13="low13_19"),
  list(tag="age >=65 (2019)", dat=b[b$age13>=65, ],  l11="low11_19", l13="low13_19"),
  list(tag="50-64 (2019)", dat=b[b$band=="50-64", ], l11="low11_19", l13="low13_19"),
  list(tag="50-64 (2025)", dat=b[b$band=="50-64", ], l11="low11_25", l13="low13_25"),
  list(tag="45-49 (2019)", dat=b[b$band=="45-49", ], l11="low11_19", l13="low13_19"))

rows <- list()
for (s in specs) {
  v <- calc(s$l11, s$l13, s$dat); q <- boot(s$l11, s$l13, s$dat)
  cat(sprintf("%-15s n=%4d low11=%3d(%4.1f%%) agree=%4.1f%% k=%.2f(%.2f-%.2f) PA=%4.1f%%(%s) NA=%4.1f%% PABAK=%.2f rev=%4.1f%%(%s) inc=%4.1f%%\n",
      s$tag, v["n"], v["LL"]+v["LN"], 100*v["prev11"], 100*v["po"], v["kappa"], q[1,1], q[2,1],
      100*v["pa"], ci(q[,3]), 100*v["na"], v["pabak"], 100*v["rev"], ci(q[,2]), 100*v["inc"]))
  rows[[length(rows)+1]] <- data.frame(
    stratum=s$tag, n=as.integer(v["n"]), low_2011=sprintf("%d (%.1f)", v["LL"]+v["LN"], 100*v["prev11"]),
    agreement=lab(v["po"]), kappa=sprintf("%.2f (%.2f-%.2f)", v["kappa"], q[1,1], q[2,1]),
    PA=sprintf("%s (%s)", lab(v["pa"]), ci(q[,3])), NAagree=lab(v["na"]),
    PABAK=sprintf("%.2f", v["pabak"]), reversion=sprintf("%s (%s)", lab(v["rev"]), ci(q[,2])),
    incidence=lab(v["inc"]), stringsAsFactors=FALSE)
}
res <- bind_rows(rows)
write.csv(res, "out/95_age_strata.csv", row.names=FALSE)
writeLines(capture.output(print(res, row.names=FALSE)), "out/95_age_strata.txt")
