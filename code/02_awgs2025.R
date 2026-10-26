# 02_awgs2025.R -- Analysis (1): AWGS 2025 age-specific cut-offs (reviewer M2)
# Compares the submitted 2019 cut-offs (male <28 / female <18 at all ages) with the
# 2025 age-specific cut-offs (50-64 y: male <34 / female <20; >=65 y: 28 / 18).
source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R")
m <- readRDS("out/jamda_full.rds")
sink("out/S1_awgs2025_report.txt", split=TRUE)
cat("==================================================================\n")
cat("SUPPLEMENTARY TABLE S1 (generated", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
cat("Analysis 1 (M2): AWGS 2025 age-specific cut-offs\n")
cat("==================================================================\n\n")

## ---------- 0. reclassification ----------
cat("--- 0. Who changes status when the cut-off changes ---\n")
r11 <- table(factor(m$low11_19, 0:1), factor(m$low11_25, 0:1), useNA="no")
r13 <- table(factor(m$low13_19, 0:1), factor(m$low13_25, 0:1), useNA="no")
cat("2011: low by 2019 =", sum(m$low11_19==1, na.rm=TRUE),
    "-> low by 2025 =", sum(m$low11_25==1, na.rm=TRUE),
    "(newly labelled low:", sum(m$low11_19==0 & m$low11_25==1, na.rm=TRUE),
    "; re-labelled normal:", sum(m$low11_19==1 & m$low11_25==0, na.rm=TRUE), ")\n")
cat("2013: low by 2019 =", sum(m$low13_19==1, na.rm=TRUE),
    "-> low by 2025 =", sum(m$low13_25==1, na.rm=TRUE),
    "(newly labelled low:", sum(m$low13_19==0 & m$low13_25==1, na.rm=TRUE),
    "; re-labelled normal:", sum(m$low13_19==1 & m$low13_25==0, na.rm=TRUE), ")\n")
# who are the newly labelled: age distribution
new <- m[m$low13_19==0 & m$low13_25==1 & !is.na(m$low13_19), ]
cat("newly labelled low at 2013: n =", nrow(new), ", age mean",
    sprintf("%.1f", mean(new$age13, na.rm=TRUE)),
    ", 65+ share", sprintf("%.1f%%", 100*mean(new$age13>=65, na.rm=TRUE)), "\n\n")

## ---------- 1. transition / stability ----------
for (v in c("19","25")) {
  t <- trans_table(m, paste0("low11_",v), paste0("low13_",v))
  tt <- t$tab
  cat(sprintf("--- 1.%s 2x2 (2011 low -> 2013 low), cut-offs %s ---\n",
              ifelse(v=="19","a","b"), ifelse(v=="19","2019","2025")))
  cat(sprintf("  never-low %d | remitted %d | incident %d | persistent %d | total %d\n",
              tt[1,1], tt[2,1], tt[1,2], tt[2,2], sum(tt)))
  cat(sprintf("  observed agreement %.3f, expected %.3f, kappa %.3f | 2011 low reverted %.1f%%\n\n",
              t$po, t$pe, t$kappa, 100*tt[1,2]/sum(tt[,2])))
}

## ---------- 2. Table 1 by group ----------
tab1 <- list()
for (v in c("19","25")) {
  b <- get_base(m, paste0("grp",v))
  for (g in GRP) {
    s <- b[b$grp==g, ]
    tab1[[length(tab1)+1]] <- data.frame(cutoff=v, group=g, n=nrow(s),
      age=sprintf("%.1f (%.1f)", mean(s$age13,na.rm=TRUE), sd(s$age13,na.rm=TRUE)),
      male=sprintf("%.1f", 100*mean(s$male,na.rm=TRUE)),
      grip13=sprintf("%.1f (%.1f)", mean(s$g13,na.rm=TRUE), sd(s$g13,na.rm=TRUE)),
      bmi=sprintf("%.1f (%.1f)", mean(s$bmi,na.rm=TRUE), sd(s$bmi,na.rm=TRUE)),
      chronic=sprintf("%.1f (%.1f)", mean(s$nchronic,na.rm=TRUE), sd(s$nchronic,na.rm=TRUE)),
      cesd=sprintf("%.1f (%.1f)", mean(s$cesd,na.rm=TRUE), sd(s$cesd,na.rm=TRUE)),
      smoke=sprintf("%.1f", 100*mean(s$smoke,na.rm=TRUE)),
      edu_high=sprintf("%.1f", 100*mean(s$edu>=4,na.rm=TRUE)), stringsAsFactors=FALSE)
  }
}
tab1 <- do.call(rbind, tab1)
cat("--- 2. Table 1 by cut-off definition ---\n")
print(tab1, row.names=FALSE)
cat("\n")

## ---------- 3. models ----------
res_all <- list(); rrr_all <- list()
for (v in c("19","25")) {
  b <- get_base(m, paste0("grp",v))
  cat(sprintf("--- 3.%s cohort (cut-offs %s): n = %d, deaths = %d, disabled = %d ---\n",
              ifelse(v=="19","a","b"), ifelse(v=="19","2019","2025"),
              nrow(b), sum(b$oc==2, na.rm=TRUE), sum(b$oc==1, na.rm=TRUE)))
  r1 <- fit_mn(b, c("age10","male"))
  r2 <- fit_mn(b, COV)
  res_all[[paste0("crude_",v)]] <- r1
  res_all[[paste0("full_",v)]  ] <- r2
  rr <- rbind(rrr_rows(r1, paste0("age+sex, cut-offs ", v)),
              rrr_rows(r2, paste0("fully adjusted, cut-offs ", v)))
  rrr_all[[v]] <- rr
  cat("  age+sex (n=", r1$n, "):\n", sep=""); print(rr[1:6, c("outcome","group","RRR","lo","hi","p")], row.names=FALSE)
  cat("  fully adjusted (n=", r2$n, "):\n", sep=""); print(rr[7:12, c("outcome","group","RRR","lo","hi","p")], row.names=FALSE)
  cat("\n")
}
rrr <- do.call(rbind, rrr_all); rownames(rrr) <- NULL
write.csv(rrr, "out/S1_rrr_19_vs_25.csv", row.names=FALSE)
write.csv(tab1, "out/S1_table1_19_vs_25.csv", row.names=FALSE)

## ---------- 4. absolute risks ----------
pr <- list()
for (v in c("19","25")) {
  p <- pred_risk(res_all[[paste0("full_",v)]])
  p$cutoff <- v; pr[[v]] <- p
}
pr <- do.call(rbind, pr); rownames(pr) <- NULL
cat("--- 4. predicted 5-year absolute risk (%), fully adjusted, 50% male ---\n")
for (v in c("19","25")) {
  cat("cut-offs", v, "\n")
  w <- reshape(pr[pr$cutoff==v, c("age","group","disabled","dead")],
               idvar="age", timevar="group", direction="wide")
  print(w, row.names=FALSE, digits=1)
}
write.csv(pr, "out/S1_absolute_risk_19_vs_25.csv", row.names=FALSE)

## ---------- 5. reversion to normal by 2015 ----------
cat("\n--- 5. persistence check (reversion by 2015, descriptive) ---\n")
for (v in c("19","25")) {
  p <- m[m[[paste0("grp",v)]] %in% "persistent" & !is.na(m$g15) & !is.na(m$male), ]
  cut <- if (v=="19") ifelse(p$male==1,28,18) else ifelse(p$male==1, ifelse(2015-p$birth<65,34,28), ifelse(2015-p$birth<65,20,18))
  cat(sprintf("  cut-offs %s: persistent with 2015 grip n = %d, normal by 2015 = %d (%.1f%%)\n",
              v, nrow(p), sum(p$g15>=cut), 100*mean(p$g15>=cut)))
}
sink()
cat("\nS1 done -> out/S1_*.csv, out/S1_awgs2025_report.txt\n")
