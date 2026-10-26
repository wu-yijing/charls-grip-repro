# 00_make_synthetic_raw.R -- generate a fully synthetic CHARLS-like raw file set
#
# Purpose: let the whole 26-script pipeline run end to end WITHOUT any CHARLS data.
# Nothing here is derived from a real respondent: every value is drawn from a
# parametric model. Distributions are calibrated only to published marginals
# (age range, sex ratio, grip cut-offs, ADL and death prevalences).
#
# Output: data/raw_synthetic/{2011,2013,2015,2018}/*.dta
# Run from the repository root:  Rscript code/00_make_synthetic_raw.R

suppressMessages({library(haven)})
set.seed(20261024)
dir.create("data/raw_synthetic", showWarnings = FALSE, recursive = TRUE)
for (y in c("2011","2013","2015","2018")) dir.create(file.path("data/raw_synthetic", y), showWarnings = FALSE)

N <- 17708
## CHARLS ID convention: the 2011 files carry an 11-character ID and conv11()
## restores the 12-character form by inserting "0" at position 10. Household in
## the first 9 digits, person index in the last 2.
hh   <- (seq_len(N) - 1) %/% 99 + 1
pp   <- (seq_len(N) - 1) %% 99 + 1
id11 <- paste0(sprintf("%09.0f", 100000000 + hh), sprintf("%02d", pp))
id12 <- paste0(substr(id11, 1, 9), "0", substr(id11, 10, 11))

## ---- sampling frame -------------------------------------------------------
ncom <- 450
com  <- sample(seq_len(ncom), N, replace = TRUE)
prov <- c("Anhui","Beijing","Chongqing","Fujian","Gansu","Guangdong","Guangxi","Guizhou",
          "Hainan","Hebei","Heilongjiang","Henan","Hubei","Hunan","Jiangsu","Jiangxi",
          "Jilin","Liaoning","Neimenggu","Ningxia","Qinghai","Shaanxi","Shandong",
          "Shanghai","Shanxi","Sichuan","Tianjin","Yunnan","Zhejiang")
com_prov <- sample(prov, ncom, replace = TRUE)
com_urb  <- rbinom(ncom, 1, 0.42)
communityID <- sprintf("%012.0f", 200000000000 + seq_len(ncom))

## ---- person-level latent structure ---------------------------------------
age11 <- pmin(85, pmax(45, round(rnorm(N, 58, 10))))
male  <- rbinom(N, 1, 0.475)
birth <- 2011 - age11

## latent grip ability (kg) and per-wave observed value
mu_g <- 32.0 + 9.0 * male - 0.22 * (age11 - 45) - 0.06 * (2013 - 2011)
A    <- rnorm(N, mu_g, 6.5)                      # stable ability
g11t <- A + rnorm(N, 0, 5.2)                     # 2013 not yet observed
g13t <- A + rnorm(N, -0.35, 5.2)                 # slight decline
g15t <- A + rnorm(N, -0.60, 5.2)

trials <- function(obs, miss_prob) {
  M <- sapply(1:4, function(k) ifelse(runif(N) < miss_prob, 993, round(obs + rnorm(N, 0, 1.0), 1)))
  M[M < 1] <- 1
  M
}

## ---- ADL / IADL (1 = no difficulty ... 4 = cannot do) ---------------------
D13 <- -3.05 + 0.050 * (age11 - 45) + 0.35 * male - 0.055 * (A - 30) + rnorm(N, 0, 0.85)
D11 <- D13 - 0.35
D18 <- D13 + 0.45

adl_items <- function(D, n_items, offset = 0) {
  out <- sapply(seq_len(n_items), function(k) {
    z <- D + offset + seq(-0.9, 0.9, length.out = n_items)[k] + rnorm(N, 0, 0.7)
    cut(z, c(-Inf, -1.6, -0.5, 0.6, Inf), labels = FALSE)
  })
  matrix(pmin(4, pmax(1, out)), nrow = N)
}

## ---- CES-D items (1-4) ----------------------------------------------------
dist <- 0.42 * (age11 - 45) / 10 + 0.55 * (D13 + 3.05) + rnorm(N, 0, 0.9)
cesd_items <- sapply(1:10, function(k) {
  z <- dist + seq(-0.8, 1.0, length.out = 10)[k] + rnorm(N, 0, 0.6)
  cut(z, c(-Inf, -0.7, 0.3, 1.3, Inf), labels = FALSE)
})

## ---- chronic conditions, smoking, education, weight -----------------------
chronic <- sapply(1:14, function(k) {
  p <- plogis(-2.0 + 0.030 * (age11 - 45) + seq(-0.9, 0.9, length.out = 14)[k])
  ifelse(runif(N) < p, 1, 2)
})
smoke   <- ifelse(runif(N) < 0.28, 1, ifelse(runif(N) < 0.15, 3, 2))
edu_raw <- sample(1:6, N, TRUE, prob = c(.30, .28, .18, .10, .08, .06))
height  <- round(pmin(200, pmax(120, rnorm(N, 162 + 9 * male, 7))), 1)
bmi_t   <- rnorm(N, 23.6 + 2.1 * (age11 - 45) / 30, 3.4)
weight  <- round(pmin(200, pmax(25, bmi_t * (height / 100)^2)), 1)

## ---- death by 2018 --------------------------------------------------------
lp   <- -3.35 + 0.055 * (age11 - 45) + 0.45 * (D13 > -0.5) + 0.25 * male
died <- rbinom(N, 1, plogis(lp))

## ---- missingness ----------------------------------------------------------
miss_adl13 <- runif(N) < 0.285          # whole block missing
miss_adl18 <- runif(N) < 0.285
miss_cesd  <- runif(N) < 0.055
miss_bmi   <- runif(N) < 0.050
miss_smoke <- runif(N) < 0.030
miss_edu   <- runif(N) < 0.020

in11 <- runif(N) < 0.790                # measured in 2011
in13 <- runif(N) < 0.745
in15 <- runif(N) < 0.620

## ---- write the 11 source files -------------------------------------------
w <- function(df, year, name) write_dta(df, file.path("data/raw_synthetic", year, paste0(name, ".dta")))

M11 <- trials(g11t, 0.02); M13 <- trials(g13t, 0.02); M15 <- trials(g15t, 0.02)
biomark <- function(M, keep, id) {
  d <- data.frame(ID = id[keep], qc003 = M[keep, 1], qc004 = M[keep, 2],
                  qc005 = M[keep, 3], qc006 = M[keep, 4], stringsAsFactors = FALSE)
  d
}
w(biomark(M11, in11, id11), "2011", "biomarkers")
w(cbind(biomark(M13, in13, id12), data.frame(
    qi002 = ifelse(miss_bmi, NA, height)[in13], ql002 = ifelse(miss_bmi, NA, weight)[in13])), "2013", "Biomarker")
w(biomark(M15, in15, id12), "2015", "Biomarker")

A11 <- adl_items(D11, 6); A13 <- adl_items(D13, 6)
I13 <- adl_items(D13, 5, offset = -0.4); A18 <- adl_items(D18, 6); I18 <- adl_items(D18, 5, offset = -0.4)
A11[runif(N) < 0.12, ] <- NA
A13[miss_adl13, ] <- NA
A18[miss_adl18, ] <- NA

h11 <- data.frame(ID = id11, as.data.frame(A11), chronic, da059 = smoke, stringsAsFactors = FALSE)
names(h11)[2:7] <- paste0("db0", 10:15)
names(h11)[8:21] <- paste0("da007_", 1:14, "_")
w(h11, "2011", "health_status_and_functioning")

h13 <- data.frame(ID = id12, as.data.frame(A13), cesd_items, stringsAsFactors = FALSE)
names(h13)[2:7] <- paste0("db0", 10:15)
names(h13)[8:17] <- sprintf("dc%03d", 9:18)
if (any(miss_cesd)) for (k in 8:17) h13[[k]][miss_cesd] <- NA
w(h13, "2013", "Health_Status_and_Functioning")

h18 <- data.frame(ID = id12, as.data.frame(A18), as.data.frame(I18), stringsAsFactors = FALSE)
names(h18)[2:7] <- paste0("db0", 10:15)
names(h18)[8:12] <- paste0("db0", 16:20)
w(h18, "2018", "Health_Status_and_Functioning")

demo <- data.frame(ID = id11, rgender = as.character(ifelse(male == 1, 1, 2)),
                   ba002_1 = birth, bd001 = ifelse(miss_edu, NA, edu_raw), stringsAsFactors = FALSE)
w(demo, "2011", "demographic_background")

w(data.frame(ID = id12, died = died, stringsAsFactors = FALSE), "2018", "Sample_Infor")

wt <- function(mu) round(pmin(6, pmax(0.05, exp(rnorm(N, log(mu), 0.45)))), 4)
w(data.frame(ID = id12, communityID = communityID[com],
             INDV_weight_ad2 = ifelse(runif(N) < 0.04, NA, wt(1)),
             INDV_weight_ad1 = ifelse(runif(N) < 0.04, NA, wt(1)),
             INDV_L_weight  = ifelse(runif(N) < 0.06, NA, wt(1)), stringsAsFactors = FALSE),
  "2013", "Weights")
w(data.frame(ID = id11, ind_weight_ad2 = ifelse(runif(N) < 0.05, NA, wt(1)), stringsAsFactors = FALSE),
  "2011", "weight")
w(data.frame(ID = id12, INDV_weight_ad2 = ifelse(runif(N) < 0.05, NA, wt(1)), stringsAsFactors = FALSE),
  "2018", "Weights")
w(data.frame(communityID = communityID, province_eng = com_prov, urban_nbs = com_urb, stringsAsFactors = FALSE),
  "2011", "PSU")

cat("synthetic raw files written to data/raw_synthetic/\n")
cat("persons:", N, "| communities:", ncom, "\n")
cat("low grip 2011 (AWGS 2019, by sex):",
    round(100 * mean(ifelse(male == 1, g11t < 28, g11t < 18)), 1), "%\n")
