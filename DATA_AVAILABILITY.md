# Data availability — what is shared and what is not

## Statement in the manuscript

> CHARLS data are available to registered users from the China Health and Retirement
> Longitudinal Study (http://charls.pku.edu.cn) under a data-use agreement that prohibits
> redistribution; the individual-level files therefore cannot be deposited or shared by the
> authors. The analysis code, a codebook, the fitted model objects behind the standardised
> risks and a deterministic script that generates a synthetic cohort of identical structure
> and missingness (seed 20261024) are openly available at Zenodo (DOI
> 10.5281/zenodo.23732082, version 1.0.4; the concept DOI 10.5281/zenodo.23677282 resolves
> to the latest version); the analysis specification is provided as Supplementary Methods.
> Estimates underlying every table and figure are reported in the manuscript and
> supplementary tables.

The quoted passage is the statement as it appears in the manuscript. The code is also browsable at
<https://github.com/wu-yijing/charls-grip-repro>, but the Zenodo record — not the repository — is the
archival copy, and the one to cite.

## Why the individual-level files cannot be posted

The CHARLS terms of agreement state, in the version every registered user accepts:

1. *Must only* use CHARLS data for academic study; any use for commercial purpose is prohibited.
2. *Must not* publish original micro data in any form; all data will be used to generate statistics
   or analytical results that can not be linked to any individual respondent.
3. *Must not* disclose, release, sell, rent, lease, loan, or otherwise grant access to the data
   covered to any person.

The Chinese text is equally explicit: 「不得将CHARLS原始微观数据以任何形式出版」、
「不得将CHARLS数据发布、转让、出售、出租、租借、借贷或以其他方式赠予他人」.

Clause 3 covers *the data covered*, not only the raw `.dta` files. An analysis file with one row per
respondent — grip values, ADL items, death status, community identifier, sampling weight — is still
covered, whether it is public or in a private repository, and whether it is posted or sent on
request. For that reason the manuscript does not offer "data available from the corresponding
author"; requests are directed to CHARLS.

## What is shared instead

| Item | Path | Why it is allowed |
|---|---|---|
| Analysis code: 40 scripts, generator, helper library | `code/` | no respondent values |
| Codebook and variable construction | `codebook/` | describes variables, not records |
| Synthetic dataset, 17,708 simulated respondents | generator only: `code/00_make_synthetic_raw.R` (the generated row-level file is deliberately not distributed) | no real person behind any row, but a per-person file is structurally indistinguishable from micro data |
| Fitted model objects (coefficients, variance) | `model/fit_store.rds` | aggregate parameters, ≈24 KB per cut-off |
| Input checksums (for version matching) | `manifest/inputs.sha256` | hashes only |
| Published estimates | manuscript and supplementary tables | aggregate |

### Why the simulated cohort is not distributed as a file

The generator ships; the CSV it produces does not. Every row is simulated, but a one-row-per-person
file carrying a respondent ID, a community identifier, a province and a sampling weight is
structurally indistinguishable from the micro data the CHARLS terms prohibit redistributing: a reader
cannot tell from the file which it is, and a journal cannot police the difference. The repository
enforces this at commit time with `tools/check_no_microdata.py`, which rejects such a file whether or
not it is synthetic. Version 1.0.0 of this package said the archival snapshot included
`jamda_synthetic.csv`. It did not. Because the file can be regenerated in one command, the wrong
statement has been removed rather than made true by shipping the file.

## What a reader can verify

1. **The pipeline runs.** `bash run_all.sh` regenerates every table and figure from the synthetic data.
2. **The numbers are the published ones.** All estimates are tabulated in the manuscript and
   supplementary tables; nothing in the paper depends on an unreported computation.
3. **The model is recoverable.** `model/fit_store.rds` holds the coefficients and variance estimates,
   so Table 3 can be recomputed without respondent data.
4. **The inputs were the standard release.** `manifest/inputs.sha256` lists the SHA-256 and row/column
   counts of each CHARLS file used; a registered user can confirm they hold the same release.

## Licensing

Code: MIT (see `LICENSE`). Synthetic data: CC0. Neither licence can override the CHARLS terms, which
continue to govern any use of the real data.
