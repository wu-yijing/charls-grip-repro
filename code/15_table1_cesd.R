source(if (file.exists("code/lib.R")) "code/lib.R" else "lib.R"); suppressMessages(library(dplyr))
m <- readRDS("out/jamda_full.rds")
b <- m %>% filter(a13==0, !is.na(grp19), !is.na(age13), !is.na(male))
for (g in c("never","remitted","incident","persistent")) {
  s <- b[b$grp19==g, ]
  cat(sprintf("%-11s %.1f (%.1f)\n", g, mean(s$cesd, na.rm=TRUE), sd(s$cesd, na.rm=TRUE)))
}
cat(sprintf("overall      %.1f (%.1f)\n", mean(b$cesd,na.rm=TRUE), sd(b$cesd,na.rm=TRUE)))
cat("missing cesd in analysis set:", sum(is.na(b$cesd)), "\n")
