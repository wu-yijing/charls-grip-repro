suppressMessages(library(haven))
rows <- list()
grab <- function(path, year, vars, domain){
  d <- read_dta(path, n_max=200)
  for (v in vars) {
    if (!(v %in% names(d))) next
    lab <- attr(d[[v]], "label"); vl <- attr(d[[v]], "labels")
    vls <- if (is.null(vl)) "" else paste(sprintf("%s=%s", vl, names(vl)), collapse="; ")
    rows[[length(rows)+1]] <<- data.frame(year=year, domain=domain, variable=v,
      label=ifelse(is.null(lab), "", lab), value_labels=vls, stringsAsFactors=FALSE)
  }
}
grab("data/raw_synthetic/2011/biomarkers.dta", 2011, c("qc003","qc004","qc005","qc006"), "grip strength")
grab("data/raw_synthetic/2013/Biomarker.dta", 2013, c("qc003","qc004","qc005","qc006"), "grip strength")
grab("data/raw_synthetic/2015/Biomarker.dta", 2015, c("qc003","qc004","qc005","qc006"), "grip strength")
grab("data/raw_synthetic/2011/health_status_and_functioning.dta", 2011, paste0("da007_",1:14,"_"), "chronic conditions")
grab("data/raw_synthetic/2011/health_status_and_functioning.dta", 2011, "da059", "smoking")
grab("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", 2013, c("db010","db011","db012","db013","db014","db015"), "ADL")
grab("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", 2013, c("db016","db017","db018","db019","db020"), "IADL")
grab("data/raw_synthetic/2018/Health_Status_and_Functioning.dta", 2018, c("db010","db011","db012","db013","db014","db015"), "ADL")
grab("data/raw_synthetic/2018/Health_Status_and_Functioning.dta", 2018, c("db016","db017","db018","db019","db020"), "IADL")
grab("data/raw_synthetic/2013/Health_Status_and_Functioning.dta", 2013, paste0("dc01",1:8), "CES-D")
out <- do.call(rbind, rows)
write.csv(out, "out/charls_codebook_check.csv", row.names=FALSE, fileEncoding="UTF-8")
cat("rows:", nrow(out), "\n")
print(table(out$domain))
