# Scripts that need the CHARLS release

Four scripts in the package read the CHARLS microdata directly, not the synthetic cohort:
`80_rev23_review.R` (hand dominance and the raw grip trials), `_probe_recon.R`
(the chair-stand item, for the reconciliation in Supplementary Table S15), and, added in
v1.0.8, `98_full_awgs_icc.R` (Supplementary Tables S37 and S38: the full AWGS 2019
possible-sarcopenia definition in all 6,052 participants and the short-interval ICC
0.80-0.99 sensitivity) and `99_full_awgs_cc.R` (the chair-stand complete-case check).
They are kept here so that the provenance of every exhibit is on record, and they are
**not** in `run_all.sh`, because no synthetic substitute exists for the variables they need.

The `98_`/`99_` scripts run from the analysis directory that contains
`out/jamda_full.rds`, the CHARLS wave files under `_x/`, and `lib.R`; they write their
summary CSVs to `out/`. Machine-specific absolute paths are deliberately absent. The
four CSVs shipped under `manifest/s37_s38/` are the archived outputs of those runs.

To run the older two scripts, place the CHARLS waves under
`data/charls/{2011,2013,2015,2018}/` in the layout `01_build.R` expects, build
`out/jamda_full.rds` from the same release, and run from the repository root with
`Rscript code/charls_only/<script>`.

Nothing in this folder contains or redistributes CHARLS data, in any form.
