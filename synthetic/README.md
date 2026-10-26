# Synthetic data

This directory is intentionally empty in version control.

`code/00_make_synthetic_raw.R` is deterministic (seed 20261024) and writes 11 `.dta` files with the
same structure, variable names, ID conventions and missingness patterns as the 2011-2018 CHARLS
releases; `code/01_build.R` then builds `out/jamda_full.rds` from them.

```bash
Rscript code/00_make_synthetic_raw.R
Rscript code/01_build.R
```

Shipping a generator rather than a 3.8 MB CSV keeps the repository honest: the file can always be
regenerated, and no one can mistake a committed table for respondent data.

**The synthetic cohort does not reproduce the published estimates.** It exists to prove that the
code runs and that every exhibit can be rebuilt; the published numbers are in the manuscript and
supplementary tables, and `model/fit_store.rds` reproduces Table 3 from the fitted parameters.

**Where the file actually lives.** The git history ships the generator only. The Zenodo snapshot
taken from a release includes `jamda_synthetic.csv` (≈3.8 MB) so that a reader without R can still
obtain the simulated cohort; that archival copy is what the manuscript cites.
