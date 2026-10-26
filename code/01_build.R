# 01_build.R -- rebuild the JAMDA analysis cohort in R
# Replicates the Python pipeline (analysis_b1.py + analysis_jamda.py) exactly.
# Gate: must reproduce n = 6,052; groups 4,943 / 385 / 462 / 262; deaths 286.
suppressMessages({library(haven); library(dplyr)})
num <- function(x){ x <- zap_labels(x); if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x) }
conv11 <- function(x) paste0(substr(x,1,9), "0", substr(x,10,11))

## ---------- grip ----------
grip <- function(path, conv=FALSE){
  d <- read_dta(path, col_select = c(ID, qc003, qc004, qc005, qc006))
  m <- sapply(d[c("qc003","qc004","qc005","qc006")], num)
  m[m %in% c(993,999)] <- NA
  m[m < 1 | m > 80] <- NA
  g <- apply(m, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm=TRUE))
  id <- as.character(d$ID); if (conv) id <- conv11(id)
  data.frame(ID = id, grip = g, stringsAsFactors = FALSE)
}
g11 <- grip("data/raw_synthetic/2011/biomarkers.dta", TRUE);  names(g11)[2] <- "g11"
g13 <- grip("data/raw_synthetic/2013/Biomarker.dta");         names(g13)[2] <- "g13"
g15 <- grip("data/raw_synthetic/2015/Biomarker.dta");         names(g15)[2] <- "g15"

## ---------- ADL / IADL ----------
adlset <- function(path, cols, thresh=3){
  d <- read_dta(path, col_select = all_of(c("ID", cols)))
  M <- sapply(d[cols], num)
  y <- as.numeric(apply(M >= thresh, 1, any, na.rm=TRUE))
  y[apply(is.na(M), 1, all)] <- NA
  data.frame(ID = as.character(d$ID), y = y, stringsAsFactors = FALSE)
}
ADL  <- c("db010","db011","db012","db013","db014","db015")
IADL <- c("db016","db017","db018","db019","db020")

a11 <- adlset("data/raw_synthetic/2011/health_status_and_functioning.dta", ADL); names(a11)[2] <- "a11"
a11$ID <- conv11(a11$ID)
a13 <- adlset("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", ADL); names(a13)[2] <- "a13"
h18 <- "data/raw_synthetic/2018/Health_Status_and_Functioning.dta"
aD  <- adlset(h18, ADL, 3);       names(aD)[2]  <- "aD"
aD4 <- adlset(h18, ADL, 4);       names(aD4)[2] <- "aD4"
aC  <- adlset(h18, c(ADL,IADL), 3); names(aC)[2] <- "aC"

## ---------- demographics, BMI, covariates ----------
demo <- read_dta("data/raw_synthetic/2011/demographic_background.dta", col_select = c(ID, rgender, ba002_1, bd001))
demo <- data.frame(ID = conv11(as.character(demo$ID)),
                   male = as.numeric(sub("^(\\d).*", "\\1", as.character(demo$rgender))) %>% {ifelse(.==1,1,ifelse(.==2,0,NA))},
                   birth = num(demo$ba002_1),
                   edu_raw = num(demo$bd001), stringsAsFactors = FALSE) %>% distinct(ID, .keep_all=TRUE)
demo$edu <- with(demo, ifelse(edu_raw==1, 0,
                       ifelse(edu_raw %in% c(2,3), 1,
                       ifelse(edu_raw==4, 2,
                       ifelse(edu_raw==5, 3, ifelse(edu_raw>=6, 4, NA))))))

b13 <- read_dta("data/raw_synthetic/2013/Biomarker.dta", col_select = c(ID, qi002, ql002))
b13 <- data.frame(ID = as.character(b13$ID), height = num(b13$qi002), weight = num(b13$ql002))
b13$height[!(b13$height>=120 & b13$height<=200)] <- NA
b13$weight[!(b13$weight>=25  & b13$weight<=200)]  <- NA
b13$bmi <- b13$weight/(b13$height/100)^2
b13$bmi[b13$bmi<12 | b13$bmi>50] <- NA
b13 <- b13 %>% distinct(ID, .keep_all=TRUE) %>% select(ID, bmi)

CHRONIC <- paste0("da007_", 1:14, "_")
h11 <- read_dta("data/raw_synthetic/2011/health_status_and_functioning.dta", col_select = all_of(c("ID","da059", CHRONIC)))
h11 <- data.frame(ID = conv11(as.character(h11$ID)), h11[CHRONIC], smoke_raw = as.character(h11$da059))
ch <- sapply(h11[CHRONIC], function(s) as.numeric(grepl("^1", as.character(num(s)))))
h11$nchronic <- rowSums(ch, na.rm=TRUE); h11$nchronic[rowSums(!is.na(ch)) < 10] <- NA
h11$smoke <- ifelse(grepl("^1", h11$smoke_raw), 1, ifelse(grepl("^2", h11$smoke_raw), 0, NA))
h11 <- h11 %>% distinct(ID, .keep_all=TRUE) %>% select(ID, nchronic, smoke)

# Full 10-item CES-D administered in 2013 (dc009-dc018), not the 8-item subset.
CESD <- c("dc009","dc010","dc011","dc012","dc013","dc014","dc015","dc016","dc017","dc018")
h13 <- read_dta("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", col_select = all_of(c("ID", CESD)))
Cm <- sapply(h13[CESD], num)
# dc013 ("hopeful about the future") and dc016 ("happy") are positively worded and must be
# reverse-scored (5 - x) so that higher scores mean more symptoms.
Cm[, c("dc013","dc016")] <- 5 - Cm[, c("dc013","dc016")]
# Prorated score on the 10-40 scale: mean of available items x 10, requiring >=8 of 10 items,
# so that partial responders are not systematically scored as less depressed.
nitem <- rowSums(!is.na(Cm))
h13 <- data.frame(ID = as.character(h13$ID), cesd = rowMeans(Cm, na.rm = TRUE) * 10)
h13$cesd[nitem < 8] <- NA
h13 <- h13 %>% distinct(ID, .keep_all=TRUE)

d18 <- read_dta("data/raw_synthetic/2018/Sample_Infor.dta", col_select = c(ID, died))
d18 <- data.frame(ID = as.character(d18$ID),
                  dead18 = as.integer(grepl("^1", as.character(num(d18$died))))) %>%
       distinct(ID, .keep_all=TRUE)

## ---------- merge ----------
m <- g11 %>% full_join(g13, by="ID") %>% full_join(g15, by="ID") %>%
  left_join(a11, by="ID") %>% left_join(a13, by="ID") %>%
  left_join(demo %>% select(ID, male, birth, edu), by="ID") %>%
  left_join(b13, by="ID") %>% left_join(h11, by="ID") %>% left_join(h13, by="ID") %>%
  left_join(aD, by="ID") %>% left_join(aD4, by="ID") %>% left_join(aC, by="ID") %>%
  left_join(d18, by="ID")
m$dead18[is.na(m$dead18)] <- 0L
m$age13 <- 2013 - m$birth
m$age11 <- 2011 - m$birth

## ---------- weights and PSU ----------
w13 <- read_dta("data/raw_synthetic/2013/Weights.dta", col_select = c(ID, communityID, INDV_weight_ad2, INDV_weight_ad1, INDV_L_weight))
w13 <- data.frame(ID = as.character(w13$ID), communityID = as.character(w13$communityID),
                  w13 = num(w13$INDV_weight_ad2), w13_raw = num(w13$INDV_weight_ad1),
                  w13_long = num(w13$INDV_L_weight)) %>% distinct(ID, .keep_all=TRUE)
w11 <- read_dta("data/raw_synthetic/2011/weight.dta", col_select = c(ID, ind_weight_ad2))
w11 <- data.frame(ID = conv11(as.character(w11$ID)), w11 = num(w11$ind_weight_ad2)) %>% distinct(ID, .keep_all=TRUE)
w18 <- read_dta("data/raw_synthetic/2018/Weights.dta", col_select = c(ID, INDV_weight_ad2))
w18 <- data.frame(ID = as.character(w18$ID), w18 = num(w18$INDV_weight_ad2)) %>% distinct(ID, .keep_all=TRUE)
psu <- read_dta("data/raw_synthetic/2011/PSU.dta", col_select = c(communityID, province_eng, urban_nbs))
psu <- data.frame(communityID = as.character(psu$communityID),
                  province = as.character(psu$province_eng),
                  urban = num(psu$urban_nbs)) %>% distinct(communityID, .keep_all=TRUE)

m <- m %>% left_join(w13, by="ID") %>% left_join(w11, by="ID") %>% left_join(w18, by="ID") %>%
  left_join(psu, by="communityID")
m$stratum <- paste(m$province, m$urban, sep="_")

## ---------- low-grip flags ----------
# 2019 AWGS cutoffs (used in the submitted manuscript)
m$low11_19 <- ifelse(!is.na(m$g11) & !is.na(m$male), as.numeric(m$g11 < ifelse(m$male==1,28,18)), NA)
m$low13_19 <- ifelse(!is.na(m$g13) & !is.na(m$male), as.numeric(m$g13 < ifelse(m$male==1,28,18)), NA)
# AWGS 2025 age-specific cutoffs: 50-64 y male <34 / female <20; >=65 y male <28 / female <18
cut25 <- function(age, male) ifelse(male==1, ifelse(age < 65, 34, 28), ifelse(age < 65, 20, 18))
m$low11_25 <- ifelse(!is.na(m$g11) & !is.na(m$male) & !is.na(m$age11),
                     as.numeric(m$g11 < cut25(m$age11, m$male)), NA)
m$low13_25 <- ifelse(!is.na(m$g13) & !is.na(m$male) & !is.na(m$age13),
                     as.numeric(m$g13 < cut25(m$age13, m$male)), NA)

grp_of <- function(l11, l13) ifelse(is.na(l11)|is.na(l13), NA_character_,
  ifelse(l11==1 & l13==1, "persistent",
  ifelse(l11==1 & l13==0, "remitted",
  ifelse(l11==0 & l13==1, "incident", "never"))))
m$grp19 <- grp_of(m$low11_19, m$low13_19)
m$grp25 <- grp_of(m$low11_25, m$low13_25)

m$age10 <- (m$age13 - 65)/10
for (k in 1:4) m[[paste0("edu",k)]] <- ifelse(is.na(m$edu), NA, as.numeric(m$edu==k))

## ---------- outcome ----------
m$oc  <- ifelse(m$dead18==1, 2, ifelse(!is.na(m$aD),  as.numeric(m$aD==1),  NA))
m$oc4 <- ifelse(m$dead18==1, 2, ifelse(!is.na(m$aD4), as.numeric(m$aD4==1), NA))
m$ocC <- ifelse(m$dead18==1, 2, ifelse(!is.na(m$aC),  as.numeric(m$aC==1),  NA))

saveRDS(m, "out/jamda_full.rds")
write.csv(m, "out/jamda_full.csv", row.names=FALSE, na="")

## ---------- verification gate ----------
base <- m %>% filter(a13 %in% 0, !is.na(grp19), !is.na(age13), !is.na(male))
cat("== GATE 1: cohort size ==\n")
cat("2013 ADL-free with grip pattern:", nrow(base), " (target 6052)\n")
cat("2018 ADL observed:", sum(!is.na(base$aD)), "(target 4453)",
    "| died:", sum(base$dead18), "(target 286)",
    "| lost:", sum(is.na(base$oc)), "(target 1313)\n")
cat("\n== GATE 2: group sizes (2019 cutoffs) ==\n")
print(table(base$grp19)[c("never","remitted","incident","persistent")])
cat("target: never 4943 / remitted 385 / incident 462 / persistent 262\n")
cat("\n== GATE 3: weights available in analysis set ==\n")
cat("w13 non-missing:", sum(!is.na(base$w13)), "| w11:", sum(!is.na(base$w11)),
    "| w18:", sum(!is.na(base$w18)), "| community:", sum(!is.na(base$communityID)),
    "| strata:", length(unique(na.omit(base$stratum))), "\n")
