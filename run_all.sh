#!/usr/bin/env bash
# run_all.sh -- run the full pipeline in order and report pass/fail per script.
# Usage:  bash run_all.sh            (from the repository root)
# Requires: R (>= 4.2) with haven/dplyr/mice/survey, python3 with python-docx.
set -u
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON="${PYTHON:-python3}"
mkdir -p out figures logs

ORDER=(
  00_make_synthetic_raw.R
  01_build.R
  02_awgs2025.R
  03_weighted.R
  04_incremental.R
  05_kappa.R
  06_weighted_mi.R
  06b_weighted_mi_main.R
  07_figure2.R
  07b_figure2.R
  08_supp_tables.py
  09_labels.R
  10_checks.R
  11_cesd_impact.R
  12_codebook.R
  13_unweighted_refs.R
  14_outcome_defs.R
  15_table1_cesd.R
  20a_cohort_facts.R
  20b_design_mi.R
  20c_table3ci_ipcw.R
  20d_sexint.R
  20e_rep_primary.R
  21_enhanced_A.R
  22_null_sim.R
  23_null_sensitivity.R
  24_n5_mi_recheck.R
  30_rev9_incremental.R
  32_rev9_boot.R
  40_share_vs_charls.R
  44_extra.R
  48_lca.R
  50_rev11_review.R
  60_rev14_spline.R
  70_rev15_decomp.R
  72_calibration.R
  73_bounding_mi.R
  74_frailty_mi.R
  75_kappa_cluster.R
  99_export_synthetic_csv.R
)

pass=0; fail=0; failed=""
for f in "${ORDER[@]}"; do
  log="logs/${f%.*}.log"
  start=$(date +%s)
  case "$f" in
    *.R)  "$RSCRIPT" "code/$f"  >"$log" 2>&1 ;;
    *.py) "$PYTHON"  "code/$f"  >"$log" 2>&1 ;;
  esac
  rc=$?
  el=$(( $(date +%s) - start ))
  if [ $rc -eq 0 ]; then
    pass=$((pass+1)); printf "PASS  %-28s %4ss\n" "$f" "$el"
  else
    fail=$((fail+1)); failed="$failed $f"; printf "FAIL  %-28s %4ss  (see %s)\n" "$f" "$el" "$log"
  fi
done

printf 'passed: %s\nfailed: %s\n' "$pass" "$fail" > logs_summary.txt
[ -n "$failed" ] && printf 'FAILED_SCRIPTS:%s\n' "$failed" >> logs_summary.txt

echo
echo "passed: $pass  failed: $fail"
if [ "$fail" -gt 0 ]; then
  echo "failed:$failed"
  exit 1
fi
exit 0
