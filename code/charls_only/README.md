# Scripts that need the CHARLS release

Two scripts in the package read the CHARLS microdata directly, not the synthetic cohort:
`80_rev23_review.R` (hand dominance and the raw grip trials) and `_probe_recon.R` (the
chair-stand item, for the reconciliation in Supplementary Table S15). They are kept here so
that the provenance of every exhibit is on record, and they are **not** in `run_all.sh`,
because no synthetic substitute exists for the variables they need.

To run them, place the CHARLS waves under `data/charls/{2011,2013,2015,2018}/` in the layout
`01_build.R` expects, build `out/jamda_full.rds` from the same release, and run from the
repository root with `Rscript code/charls_only/<script>`.

Nothing in this folder contains or redistributes CHARLS data, in any form.
