# Synthetic data

This directory ships the generator, not its output.

`code/00_make_synthetic_raw.R` is deterministic (seed 20261024) and writes 11 `.dta` files with the
same structure, variable names, ID conventions and missingness patterns as the 2011-2018 CHARLS
releases; `code/01_build.R` then builds `out/jamda_full.rds` from them, and
`code/99_export_synthetic_csv.R` writes the same cohort out as a local CSV.

```bash
Rscript code/00_make_synthetic_raw.R
Rscript code/01_build.R
Rscript code/99_export_synthetic_csv.R     # optional, writes synthetic/jamda_synthetic.csv
```

**Why the CSV is not committed.** Every row is simulated, but a one-row-per-person file carrying an
ID, a community identifier, a province and a sampling weight is structurally indistinguishable from
the micro data the CHARLS terms of use prohibit redistributing, and `tools/check_no_microdata.py`
rejects it on exactly that ground. The guard is deliberate and is not relaxed for synthetic files.

Version 1.0.0 of this package stated that the archival snapshot included `jamda_synthetic.csv`. It
did not; the statement was wrong, and it has been removed rather than made true by shipping the file.

**The synthetic cohort does not reproduce the published estimates.** It exists to prove that the
code runs and that every exhibit can be rebuilt; the published numbers are in the manuscript and
supplementary tables, and `model/fit_store.rds` reproduces Table 3 from the fitted parameters.
