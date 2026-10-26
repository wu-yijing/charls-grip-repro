#!/usr/bin/env python3
"""Block staged files larger than 10 MiB (GitHub warns at 50 MiB, blocks at 100 MiB)."""
import os, subprocess, sys

MAX = 10 * 1024 * 1024
out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
                     capture_output=True, text=True, check=False).stdout.split()
bad = [(f, os.path.getsize(f)) for f in out if os.path.isfile(f) and os.path.getsize(f) > MAX]
if not bad:
    print("check_size: OK (%d file(s))" % len(out)); sys.exit(0)
print("\nCOMMIT BLOCKED — files exceed 10 MiB:/n")
for f, s in bad:
    print(f"  {f}  {s/1048576:.1f} MiB")
print("\nLarge binaries do not belong in a source repository; ship a script that regenerates them.")
sys.exit(1)
