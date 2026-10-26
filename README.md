# Reproducibility package

[![pipeline](https://github.com/wu-yijing/charls-grip-repro/actions/workflows/repro.yml/badge.svg)](https://github.com/wu-yijing/charls-grip-repro/actions/workflows/repro.yml)

**Reproducibility of low grip strength classification and the incremental value of a prior reading beyond the current value: a five-year cohort study of 6,052 Chinese adults**

This repository contains the 26 analysis scripts (plus a synthetic-data generator and a shared
helper library), the analysis specification, a synthetic dataset of identical structure, and the
fitted model objects for the manuscript. It contains **no CHARLS
data** — see [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for why, and for what is shared instead.

---

## Quick start

```bash
# 1. generate a synthetic raw file set (11 .dta files, ~30 MB; nothing real)
Rscript code/00_make_synthetic_raw.R

# 2. run all 27 steps in dependency order (26 analysis scripts + the generator)
bash run_all.sh
```

`run_all.sh` prints one `PASS`/`FAIL` line per step, writes `logs_summary.txt`, and exits 1 if any
step failed (so CI fails loudly). Every step writes its own log to `logs/`.

Requirements: R ≥ 4.2 (`haven`, `dplyr`, `nnet`, `sandwich`, `pROC`, `mice`, `survey`) and
Python 3 with `python-docx` for the supplementary-table builder.

---

## What the synthetic data are, and are not

`code/00_make_synthetic_raw.R` generates 17,708 simulated respondents with a latent grip ability,
correlated ADL and death outcomes, and the same file structure, variable names, ID conventions and
missingness patterns as the 2011–2018 CHARLS releases. `code/01_build.R` then builds the analysis
cohort from those files exactly as it does from the real ones.

* **It is enough to verify** that the code runs, that every table and figure can be regenerated,
  and that the analytic decisions are implemented as described.
* **It is not enough to reproduce the published numbers.** Estimates from the synthetic data differ
  from the paper by construction. The published estimates are in the manuscript and the
  supplementary tables, and the fitted model objects in `model/` reproduce the standardized risks
  in Table 3 exactly.

Synthetic cohort from the shipped seed: n = 5,180 (4,630 never low / 209 remitted / 231 incident /
110 persistent). The real cohort is n = 6,052 (4,943 / 385 / 462 / 262).

---

## Script map

| Script | What it does | Consumes | Produces |
|---|---|---|---|
| `00_make_synthetic_raw.R` | synthetic raw file set | — | `data/raw_synthetic/**` |
| `01_build.R` | build the cohort, flags, outcome | `data/raw_synthetic/**` | `out/jamda_full.rds`, `.csv` |
| `02_awgs2025.R` | 2019 vs 2025 cut-offs | `jamda_full.rds` | `out/S1_*` |
| `03_weighted.R` | weighted / design-based models | `jamda_full.rds` | `out/S2_multinom_weighted.csv`, `S2_design_binomial.csv` |
| `04_incremental.R` | AUC, NRI, IDI, decision curves | `jamda_full.rds` | `out/S3_incremental.csv` |
| `05_kappa.R` | agreement statistics | `jamda_full.rds` | console |
| `06_weighted_mi.R`, `06b_weighted_mi_main.R` | multiply imputed main models | `jamda_full.rds` | `out/S2_weighted_MI_*.csv` |
| `07_figure2.R`, `07b_figure2.R` | Figure 2 panels | `out/S2_weighted_MI_absolute_risk.csv`, `jamda_full.rds` | `figures/*.png` |
| `08_supp_tables.py` | assemble S1–S3 | `out/S*.csv` | `out/Supplementary_Tables_S1-S3.docx` |
| `09_labels.R`, `12_codebook.R` | variable labels / codebook | raw `.dta` | console |
| `10_checks.R`, `11_cesd_impact.R` | CES-D construction checks | raw + built | console, `jamda_full_fix.rds` |
| `13_unweighted_refs.R` | unweighted reference models | `jamda_full.rds` | console |
| `14_outcome_defs.R` | alternative outcome definitions, sex interaction | `jamda_full.rds` | console |
| `15_table1_cesd.R` | Table 1 CES-D row | `jamda_full.rds` | console |
| `20a_cohort_facts.R` | cohort flow, numbers quoted in text | `jamda_full.rds` | `out/20a_cohort_facts.txt` |
| `20b_design_mi.R` | design-based + MI sensitivity, model store | `jamda_full.rds` | `out/S5_*`, `out/20b_fit_store.rds` |
| `20c_table3ci_ipcw.R` | Table 3 confidence intervals | `jamda_full.rds`, `20b_fit_store.rds` | `out/S6_absolute_risk_ci.csv` |
| `20d_sexint.R` | sex × pattern interaction | `jamda_full.rds` | `out/20d_sexint.txt` |
| `20e_rep_primary.R` | replication of the primary model | `jamda_full.rds` | console |
| `21_enhanced_A.R` | κ CI, PABAK, positive agreement, SEM, MDC | `jamda_full.rds` | `out/21_enhanced_A_v2.txt` |
| `22_null_sim.R` | measurement-error-only simulation | `jamda_full.rds` | console |
| `23_null_sensitivity.R` | sensitivity of the noise share | `jamda_full.rds` | console |
| `24_n5_mi_recheck.R` | ordinal test and DCA in imputed data | `jamda_full.rds` | `out/24_n5.txt` |

`manifest/outputs.csv` maps every manuscript exhibit to the script that produces it.

---

## Data access

CHARLS data are available to registered users from the China Health and Retirement Longitudinal
Study (<http://charls.pku.edu.cn>). Registration is individual, and the data-use agreement prohibits
redistribution by users — including derived individual-level files. Requests for the data should go
to CHARLS, not to the authors.

The fitted model objects in `model/` (≈24 KB each) contain coefficients and variance estimates only;
they let a reader recompute the standardized risks in Table 3 without any respondent-level input.

---

## Guard against accidental commits

Two independent backstops block respondent-level data from entering history — a rule that cannot be
undone once pushed:

```bash
git config core.hooksPath .githooks     # no extra dependencies
# or, if you use pre-commit:
pre-commit install
```

`tools/check_no_microdata.py` rejects survey-data extensions, any CSV with identifier or design
columns (`ID`, `communityID`, `province`, `stratum`, `urban`, weights), any table over 1,000 rows,
and any file over 10 MiB. To see it work:

```bash
cp /dev/null test.dta && git add -f test.dta && git commit -m "test"   # must fail
git rm --cached test.dta
```

---

## Citation

See `CITATION.cff`. Code is MIT-licensed; the synthetic dataset is released under CC0.

## Contact

Yijing Wu — Department of Endocrinology, The First Affiliated Hospital of Anhui Medical University.
