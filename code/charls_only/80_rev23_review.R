# NOT RUNNABLE ON THE SYNTHETIC COHORT -- see charls_only/README.md.
# This script reads the 2011 and 2013 CHARLS biomarker files directly (hand dominance, qc002,
# and the raw grip trials), so it needs the CHARLS release placed under data/charls/.
# It produces the dominant-hand analysis (Supplementary Table S22 Panel A), the 2015 retention
# by pattern (Table S11 note) and the age specifications of the mortality model (S22 Panel B/C).
# 80_rev23_review.R -- analyses requested by the 2026-10-28 review
#  A. M1  low-grip definition: AWGS (better of two trials, dominant hand) vs the
#         submitted maximum of four trials; agreement indices and the prognostic contrast
#  B. M3  retention at the 2015 wave by pattern, and the non-nested "second reading" contrast
#  C. M7  functional form of age in the mortality model (linear / 5-year bands / spline / age >= 65)
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
suppressMessages({library(haven); library(dplyr); library(nnet); library(splines)})

m <- readRDS("out/jamda_full.rds")

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

num <- function(x){ x <- zap_labels(x); if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x) }
conv11 <- function(x) paste0(substr(x,1,9), "0", substr(x,10,11))
clean <- function(v){ v <- num(v); v[v %in% c(993,999)] <- NA; v[v < 1 | v > 80] <- NA; v }

## ---------- A. grip metrics from the four trials ----------
# labels: qc003 Left-1, qc004 Right-1, qc005 Left-2, qc006 Right-2
read_grip <- function(path, conv){
  d <- read_dta(path, col_select = c(ID, qc002, qc003, qc004, qc005, qc006))
  L1 <- clean(d$qc003); R1 <- clean(d$qc004); L2 <- clean(d$qc005); R2 <- clean(d$qc006)
  dom <- num(d$qc002)
  mx <- function(a, b) ifelse(is.na(a), b, ifelse(is.na(b), a, pmax(a, b)))
  L <- mx(L1, L2); R <- mx(R1, R2)
  data.frame(
    ID       = if (conv) conv11(as.character(d$ID)) else as.character(d$ID),
    g_max4   = mx(L, R),
    g_dom    = ifelse(dom == 1, R, ifelse(dom == 2, L, ifelse(dom == 3, mx(L, R), NA))),
    g_mean   = ifelse(is.na(L) | is.na(R), mx(L, R), (L + R) / 2),
    g_left   = L,
    g_right  = R,
    dom      = dom,
    stringsAsFactors = FALSE)
}
G11 <- read_grip("data/charls/2011/biomarkers.dta", TRUE)
G13 <- read_grip("data/charls/2013/Biomarker.dta",    FALSE)
names(G11)[-1] <- paste0(names(G11)[-1], "_11"); names(G13)[-1] <- paste0(names(G13)[-1], "_13")
m <- m %>% select(-any_of(c("g_max4_11","g_dom_11","g_mean_11","g_left_11","g_right_11","dom_11",
                            "g_max4_13","g_dom_13","g_mean_13","g_left_13","g_right_13","dom_13"))) %>%
  left_join(G11, by = "ID") %>% left_join(G13, by = "ID")

cat("== completeness of the four grip definitions, both waves ==\n")
for (k in c("max4","dom","mean","left","right"))
  cat(sprintf("  %-6s 2011 non-missing %5d | 2013 %5d | both %5d\n", k,
      sum(!is.na(m[[paste0("g_",k,"_11")]])), sum(!is.na(m[[paste0("g_",k,"_13")]])),
      sum(!is.na(m[[paste0("g_",k,"_11")]]) & !is.na(m[[paste0("g_",k,"_13")]]))))

agree_row <- function(l11, l13){
  ok <- !is.na(l11) & !is.na(l13)
  a <- sum(l11[ok] == 1 & l13[ok] == 1)
  b <- sum(l11[ok] == 1 & l13[ok] == 0)
  cc<- sum(l11[ok] == 0 & l13[ok] == 1)
  dd<- sum(l11[ok] == 0 & l13[ok] == 0)
  n <- a + b + cc + dd
  po <- (a + dd)/n
  p11 <- (a + b)/n; p13 <- (a + cc)/n
  pe <- p11*p13 + (1-p11)*(1-p13)
  data.frame(n = n, prev11 = 100*p11, prev13 = 100*p13, po = 100*po, pe = 100*pe,
             kappa = (po-pe)/(1-pe), pabak = 2*po - 1,
             pa = 100*2*a/(2*a + b + cc), na = 100*2*dd/(2*dd + b + cc),
             reversion = 100*b/(a+b), incidence = 100*cc/(cc+dd),
             confirm = 100*a/(a+b), stringsAsFactors = FALSE)
}

resA <- list()
for (k in c("max4","dom","mean","left","right")) {
  g11 <- m[[paste0("g_",k,"_11")]]; g13 <- m[[paste0("g_",k,"_13")]]
  l11 <- ifelse(is.na(g11) | is.na(m$male), NA, as.numeric(g11 < ifelse(m$male == 1, 28, 18)))
  l13 <- ifelse(is.na(g13) | is.na(m$male), NA, as.numeric(g13 < ifelse(m$male == 1, 28, 18)))
  keep <- m$a13 %in% 0 & !is.na(m$age13) & !is.na(m$male)
  r <- agree_row(ifelse(keep, l11, NA), ifelse(keep, l13, NA))
  r$def <- k; resA[[k]] <- r
}
resA <- do.call(rbind, resA)
cat("\n== A1. agreement of the classification under five grip definitions ==\n")
print(resA[, c("def","n","prev11","prev13","po","kappa","pabak","pa","na","reversion","incidence","confirm")],
      row.names = FALSE, digits = 3)
write.csv(resA, "out/rev23_A1_definitions.csv", row.names = FALSE)

## A2. prognostic contrast under each definition (complete case, weighted, full adjustment)
base_for <- function(k){
  g11 <- m[[paste0("g_",k,"_11")]]; g13 <- m[[paste0("g_",k,"_13")]]
  l11 <- ifelse(is.na(g11) | is.na(m$male), NA, as.numeric(g11 < ifelse(m$male == 1, 28, 18)))
  l13 <- ifelse(is.na(g13) | is.na(m$male), NA, as.numeric(g13 < ifelse(m$male == 1, 28, 18)))
  grp <- ifelse(is.na(l11)|is.na(l13), NA,
         ifelse(l11==1 & l13==1, "persistent", ifelse(l11==1 & l13==0, "remitted",
         ifelse(l11==0 & l13==1, "incident", "never"))))
  b <- m %>% mutate(grp = grp) %>% filter(a13 %in% 0, !is.na(grp), !is.na(age13), !is.na(male))
  b$grp <- relevel(factor(b$grp, levels = GRP), ref = "never")
  b$y   <- factor(b$oc, levels = 0:2, labels = c("alive_free","disabled","dead"))
  b$wsc <- b$w13 / mean(b$w13, na.rm = TRUE)
  b
}
fit_full <- function(b, covars = COV, ref = "never"){
  b$grp <- relevel(b$grp, ref = ref)
  d <- b[complete.cases(b[, c("y", covars, "grp", "wsc", "communityID")]), ]
  f <- as.formula(paste("y ~", paste(c(covars, "grp"), collapse = "+")))
  fit <- multinom(f, data = d, weights = d$wsc, trace = FALSE, maxit = 500)
  X <- cbind(1, as.matrix(d[, covars]))
  for (g in setdiff(GRP, ref)) X <- cbind(X, as.numeric(d$grp == g))
  colnames(X) <- c("const", covars, paste0("grp", setdiff(GRP, ref)))
  bb <- matrix(coef(fit), nrow = 2, byrow = FALSE)
  V <- mn_vcov(bb, X, as.integer(d$y) - 1, d$wsc, d$communityID)$V
  list(fit = fit, beta = bb, V = V, cols = colnames(X), n = nrow(d), data = d)
}
get_rrr <- function(res, g, out = 2){
  j <- which(res$cols == paste0("grp", g)); b <- res$beta[out, j]; s <- sqrt(res$V[(out-1)*length(res$cols)+j, (out-1)*length(res$cols)+j])
  data.frame(RRR = exp(b), lo = exp(b - 1.96*s), hi = exp(b + 1.96*s), p = 2*pnorm(-abs(b/s)))
}

resA2 <- list()
for (k in c("max4","dom","mean","left","right")) {
  b <- base_for(k); fit <- fit_full(b)
  r <- get_rrr(fit, "persistent")
  r$def <- k; r$n <- fit$n
  resA2[[k]] <- r
}
resA2 <- do.call(rbind, resA2)
cat("\n== A2. persistent-vs-never mortality RRR under each definition (complete case, weighted) ==\n")
print(resA2[, c("def","n","RRR","lo","hi","p")], row.names = FALSE, digits = 3)
write.csv(resA2, "out/rev23_A2_definitions_rrr.csv", row.names = FALSE)

## ---------- B. retention at 2015 and the non-nested contrast ----------
b <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
b$low15 <- ifelse(is.na(b$g15) | is.na(b$male), NA, as.numeric(b$g15 < ifelse(b$male == 1, 28, 18)))
cat("\n== B1. 2015 grip availability by 2011-2013 pattern ==\n")
tabB <- b %>% group_by(grp = factor(grp19, levels = GRP)) %>%
  summarise(n = n(), with15 = sum(!is.na(low15)), pct = 100*mean(!is.na(low15)),
            age = mean(age13, na.rm = TRUE), male = 100*mean(male, na.rm = TRUE), .groups = "drop")
print(as.data.frame(tabB), row.names = FALSE, digits = 3)
write.csv(tabB, "out/rev23_B1_retention.csv", row.names = FALSE)

cat("\n== B2. among those low in 2013 with a 2015 reading: persistence by 2011 status ==\n")
d13 <- b %>% filter(low13_19 == 1, !is.na(low15))
tb <- table(factor(d13$grp19, levels = c("persistent","incident")), d13$low15)
print(tb)
cat(sprintf("  persistent: %d/%d = %.1f%%\n", tb[1,2], sum(tb[1,]), 100*tb[1,2]/sum(tb[1,])))
cat(sprintf("  incident  : %d/%d = %.1f%%\n", tb[2,2], sum(tb[2,]), 100*tb[2,2]/sum(tb[2,])))
cat(sprintf("  marginal  : %d/%d = %.1f%%\n", sum(tb[,2]), sum(tb), 100*sum(tb[,2])/sum(tb)))
write.csv(as.data.frame(tb), "out/rev23_B2_persistence.csv")

## ---------- C. functional form of age ----------
b <- base_for("max4")
lo <- c(45,50,55,60,65,70,75,80,85,90); hi <- c(50,55,60,65,70,75,80,85,90,200)
band_cols <- character(0)
for (j in 2:length(lo)) {
  nmj <- paste0("agb", lo[j])
  b[[nmj]] <- as.numeric(b$age13 >= lo[j] & b$age13 < hi[j])
  band_cols <- c(band_cols, nmj)
}
NS  <- ns(b$age13, df = 3)
for (j in 1:3) b[[paste0("ns", j)]] <- NS[, j]
base_cov <- setdiff(COV, "age10")

cat("\n== C1. persistent-vs-never mortality RRR under four age specifications (complete case, weighted) ==\n")
specs <- list(
  "linear age10 (submitted)"    = list(data = b, cov = COV),
  "5-year age bands"            = list(data = b, cov = c(band_cols, base_cov)),
  "natural spline df=3"         = list(data = b, cov = c("ns1","ns2","ns3", base_cov)),
  "age >= 65, linear"           = list(data = subset(b, age13 >= 65), cov = COV)
)
resC <- list(); resC2 <- list()
for (nm in names(specs)) {
  s1 <- specs[[nm]]
  fit  <- fit_full(s1$data, covars = s1$cov)
  fit2 <- fit_full(s1$data, covars = s1$cov, ref = "incident")
  r  <- get_rrr(fit,  "persistent"); r$spec <- nm;  r$n <- fit$n;  r$deaths <- sum(fit$data$y == "dead")
  r2 <- get_rrr(fit2, "persistent"); r2$spec <- nm; r2$n <- fit2$n; r2$deaths <- sum(fit2$data$y == "dead")
  resC[[nm]] <- r; resC2[[nm]] <- r2
}
resC <- do.call(rbind, resC); resC2 <- do.call(rbind, resC2)
print(resC[, c("spec","n","deaths","RRR","lo","hi","p")], row.names = FALSE, digits = 3)
cat("\n== C2. persistent vs incident, mortality (incident as reference) ==\n")
print(resC2[, c("spec","n","deaths","RRR","lo","hi","p")], row.names = FALSE, digits = 3)
write.csv(rbind(cbind(panel = "C1_vs_never", resC), cbind(panel = "C2_vs_incident", resC2)),
          "out/rev23_C_age_form.csv", row.names = FALSE)

cat("\n== C3. age distribution by pattern (analytic cohort) ==\n")
print(b %>% group_by(grp) %>% summarise(n = n(), mean = mean(age13), sd = sd(age13),
      p65 = 100*mean(age13 >= 65), deaths = sum(y == "dead"), .groups = "drop") %>% as.data.frame(),
      row.names = FALSE, digits = 3)
