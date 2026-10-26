#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Guard: refuse to commit anything that could be respondent-level CHARLS data.

CHARLS terms of use prohibit publishing, disclosing or granting access to the
micro data in any form, including derived individual-level files. This hook is a
mechanical backstop for a rule that cannot be fixed after the fact: once a file
is pushed it stays in history, in forks and in caches.

Rules
  1. Forbidden extensions anywhere except the allow-listed model/ directory:
     .dta .sav .sas7bdat .xpt .rds .rdata .parquet .feather
  2. Any .csv/.tsv whose header contains an identifier or design column
     (ID, communityID, province, stratum, urban, ...) or that has > 1000 rows.
  3. Any single file larger than 10 MiB.

Usage
  python3 tools/check_no_microdata.py --staged     # pre-commit
  python3 tools/check_no_microdata.py path ...     # ad hoc
"""
import os
import subprocess
import sys

MAX_BYTES = 10 * 1024 * 1024
FORBIDDEN_EXT = {".dta", ".sav", ".sas7bdat", ".xpt", ".rds", ".rdata", ".parquet", ".feather"}
TABULAR_EXT = {".csv", ".tsv", ".txt"}
ID_COLUMNS = {
    "id", "communityid", "community", "province", "stratum", "urban", "psu",
    "hhid", "householdid", "w11", "w13", "w18", "w13_raw", "w13_long",
}
ALLOW_DIRS = (os.path.normpath("model"),)


def staged_files():
    out = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
        capture_output=True, text=True, check=False,
    )
    return [f for f in out.stdout.splitlines() if f.strip()]


def is_allowed_path(path):
    norm = os.path.normpath(path)
    return any(norm.startswith(d + os.sep) for d in ALLOW_DIRS)


def check_csv(path):
    """Return a reason string when the table looks respondent-level."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().strip().split(",")
            rows = sum(1 for _ in fh)
    except OSError as exc:
        return f"unreadable ({exc})"
    cols = {c.strip().strip('"').lower() for c in header}
    hit = sorted(cols & ID_COLUMNS)
    if hit:
        return f"identifier columns present ({', '.join(hit)})"
    if rows > 1000:
        return f"{rows:,} data rows (> 1000) — respondent-level by default"
    return None


def check(path):
    if not os.path.isfile(path):
        return None
    if is_allowed_path(path):
        return None
    ext = os.path.splitext(path)[1].lower()
    if ext in FORBIDDEN_EXT:
        return "binary survey-data format"
    size = os.path.getsize(path)
    if size > MAX_BYTES:
        return f"{size / 1048576:.1f} MiB (> 10 MiB)"
    if ext in TABULAR_EXT:
        reason = check_csv(path)
        if reason:
            return reason
    return None


def main(argv):
    files = staged_files() if "--staged" in argv else [a for a in argv if not a.startswith("-")]
    bad = [(f, check(f)) for f in files]
    bad = [(f, r) for f, r in bad if r]
    if not bad:
        print("check_no_microdata: OK (%d file(s) scanned)" % len(files))
        return 0
    print("\nCOMMIT BLOCKED — these files look like respondent-level data:\n")
    for f, r in bad:
        print(f"  {f}\n      reason: {r}")
    print(
        "\nCHARLS terms of use prohibit publishing or granting access to the micro data.\n"
        "Remove the file from the index, add it to .gitignore, and ship code, a\n"
        "codebook, the synthetic dataset or the model objects instead.\n"
        "See DATA_AVAILABILITY.md for what may be published.\n\n"
        "If a file is a false positive, move it under model/ (allow-listed) or adjust\n"
        "tools/check_no_microdata.py deliberately, and record why in the PR.\n"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
