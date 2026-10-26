suppressMessages(library(haven))
show <- function(path, vars, tag){
  d <- read_dta(path, n_max=50)
  cat("\n=====", tag, "=====\n")
  for (v in vars) {
    if (!(v %in% names(d))) { cat(sprintf("%-12s <NOT FOUND>\n", v)); next }
    lab <- attr(d[[v]], "label"); vl <- attr(d[[v]], "labels")
    cat(sprintf("%-12s label: %s\n", v, ifelse(is.null(lab), "(none)", lab)))
    if (!is.null(vl)) {
      n <- min(length(vl), 6)
      cat("             values:", paste(sprintf("%s=%s", vl[1:n], names(vl)[1:n]), collapse=" | "), "\n")
    }
  }
}
show("data/raw_synthetic/2011/biomarkers.dta", c("qc003","qc004","qc005","qc006","qi002","ql002"), "2011 biomarker (grip/height/weight)")
show("data/raw_synthetic/2011/health_status_and_functioning.dta",
     c(paste0("da007_",1:14,"_"), "da059"), "2011 HSF (14 chronic conditions, smoking)")
show("data/raw_synthetic/2013/Health_Status_and_Functioning.dta",
     c("db010","db011","db012","db013","db014","db015","db016","db017","db018","db019","db020",
       paste0("dc01",1:9), "dc020"), "2013 HSF (ADL, IADL, CES-D)")
show("data/raw_synthetic/2018/Health_Status_and_Functioning.dta",
     c("db010","db011","db012","db013","db014","db015","db016","db017","db018","db019","db020"), "2018 HSF (ADL, IADL)")
