# 24_n5_mi_recheck.R -- ordinal-constraint LR test and decision curve repeated in the imputed data
suppressMessages({library(nnet); library(dplyr); library(mice)})
m <- readRDS("out/jamda_full.rds")

base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
base$y   <- factor(base$oc, levels = 0:2, labels = c("alive_free","disabled","dead"))
base$grp <- relevel(factor(base$grp19, levels = c("never","remitted","incident","persistent")), ref = "never")
base$nr  <- base$low11_19 + base$low13_19
base$edu_f <- factor(base$edu)
base$grip5 <- base$g13 / 5

COV <- c("age10","male","bmi","nchronic","cesd","smoke")
d0 <- base[, c(COV, "edu_f", "grp", "y", "nr", "grip5", "g13", "communityID")]

set.seed(20261024)
imp <- mice(d0, m = 20, printFlag = FALSE)

cat("################ ORDINAL-CONSTRAINT LR TEST (imputed)\n")
chis <- numeric(20)
for (i in 1:20) {
  di <- complete(imp, i)
  di$y <- factor(di$y, levels = c("alive_free","disabled","dead"))
  f4 <- as.formula(paste("y ~", paste(c(COV, "edu_f", "grp"), collapse = "+")))
  f2 <- as.formula(paste("y ~", paste(c(COV, "edu_f", "nr"),  collapse = "+")))
  m4 <- multinom(f4, data = di, trace = FALSE, maxit = 500)
  m2 <- multinom(f2, data = di, trace = FALSE, maxit = 500)
  chis[i] <- 2 * (as.numeric(logLik(m4)) - as.numeric(logLik(m2)))
}
cat(sprintf("mean chi-square %.3f | range %.3f-%.3f | df 4 | P (mean) %.4f\n",
            mean(chis), min(chis), max(chis), pchisq(mean(chis), 4, lower.tail = FALSE)))

cat("\n################ DECISION CURVE (imputed)\n")
ths <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
nbA <- nbB <- matrix(NA, 20, length(ths))
nb <- function(p, y, pt) mean(p >= pt & y == 1) - mean(p >= pt & y == 0) * pt / (1 - pt)
for (i in 1:20) {
  di <- complete(imp, i)
  di <- di[di$y != "dead", ]; di$dis <- as.numeric(di$y == "disabled")
  gA <- glm(as.formula(paste("dis ~", paste(c(COV, "edu_f", "grip5"), collapse = "+"))),
            data = di, family = binomial)
  gB <- glm(as.formula(paste("dis ~", paste(c(COV, "edu_f", "grip5", "grp"), collapse = "+"))),
            data = di, family = binomial)
  pA <- fitted(gA); pB <- fitted(gB)
  for (j in seq_along(ths)) { nbA[i, j] <- nb(pA, di$dis, ths[j]); nbB[i, j] <- nb(pB, di$dis, ths[j]) }
}
res <- data.frame(threshold = ths,
                  modelA = colMeans(nbA), modelB = colMeans(nbB),
                  diff = colMeans(nbB - nbA))
print(round(res, 4), row.names = FALSE)
cat(sprintf("\nmax gain across 0.03-0.40 (imputed): %.4f\n",
            max(sapply(seq(0.03, 0.40, 0.01), function(pt) {
              d1 <- complete(imp, 1); d1 <- d1[d1$y != "dead", ]; d1$dis <- as.numeric(d1$y == "disabled")
              gA <- glm(as.formula(paste("dis ~", paste(c(COV,"edu_f","grip5"), collapse="+"))), data=d1, family=binomial)
              gB <- glm(as.formula(paste("dis ~", paste(c(COV,"edu_f","grip5","grp"), collapse="+"))), data=d1, family=binomial)
              nb(fitted(gB), d1$dis, pt) - nb(fitted(gA), d1$dis, pt)
            }))))
