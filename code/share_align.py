# Ported for the reproducibility package: run from the repository root.
# Input: out/jamda_full.csv (built by 01_build.R). Percentile of the AWGS 2019 cut-off in
# the Chinese and the European (SHARE) grip distributions, Supplementary Table S10 Panel E.
# -*- coding: utf-8 -*-
"""CHARLS vs SHARE: where does the AWGS 2019 cut-off sit in the grip distribution?"""
import csv, math

CSV = "out/jamda_full.csv"
CUT = {"men": 28.0, "women": 18.0}
BANDS = [(50,54),(55,59),(60,64),(65,69),(70,74),(75,79),(80,84),(85,200)]

# SHARE, Europe pooled, standing, absolute kg (GeroScience 2025; Tables 3 and 4)
PCTL = [5,10,20,30,40,50,60,70,80,90,95]
SHARE = {
 ("men",(50,54)): [30,35,39,41,44,46,48,51,54,58,62],
 ("men",(55,59)): [30,34,37,40,42,45,47,50,52,56,59],
 ("men",(60,64)): [27,31,36,40,42,44,46,48,51,55,57],
 ("men",(65,69)): [27,31,35,37,40,42,44,46,49,52,55],
 ("men",(70,74)): [25,28,32,35,37,39,41,44,46,49,52],
 ("men",(75,79)): [22,25,29,31,34,36,38,40,43,46,48],
 ("men",(80,84)): [20,22,27,30,31,33,35,38,40,43,45],
 ("men",(85,200)):[17,20,23,25,28,29,31,33,36,40,43],
 ("women",(50,54)): [19,21,24,26,28,29,31,32,34,37,39],
 ("women",(55,59)): [19,21,24,25,27,28,30,32,33,36,38],
 ("women",(60,64)): [19,21,24,25,26,28,29,31,32,35,37],
 ("women",(65,69)): [17,19,22,24,25,26,28,29,31,34,36],
 ("women",(70,74)): [15,17,20,22,23,25,26,28,29,32,34],
 ("women",(75,79)): [13,15,18,20,21,23,24,25,27,30,32],
 ("women",(80,84)): [12,14,17,18,20,21,22,24,25,28,30],
 ("women",(85,200)):[10,12,14,16,17,19,20,22,23,25,27],
}
SHARE_N = {("men",(50,54)):852,("men",(55,59)):2474,("men",(60,64)):3502,("men",(65,69)):3976,
           ("men",(70,74)):3826,("men",(75,79)):2626,("men",(80,84)):1574,("men",(85,200)):790,
           ("women",(50,54)):1527,("women",(55,59)):3409,("women",(60,64)):4509,("women",(65,69)):4831,
           ("women",(70,74)):4389,("women",(75,79)):3077,("women",(80,84)):1882,("women",(85,200)):970}

def pctl_pos(vals, cum, x):
    """position of x on a percentile scale, linear interpolation between tabulated points."""
    if x < vals[0]:
        return None          # below P5
    if x > vals[-1]:
        return 99.0
    for i in range(len(vals) - 1):
        if vals[i] <= x <= vals[i+1]:
            if vals[i+1] == vals[i]:
                return cum[i]
            frac = (x - vals[i]) / (vals[i+1] - vals[i])
            return cum[i] + frac * (cum[i+1] - cum[i])
    return 99.0

def band_of(age):
    for lo, hi in BANDS:
        if lo <= age <= hi:
            return (lo, hi)
    return None

rows = []
with open(CSV, encoding="utf8") as f:
    for r in csv.DictReader(f):
        try:
            g = float(r["g11"]); age = float(r["age11"]); male = int(float(r["male"]))
        except (ValueError, TypeError, KeyError):
            continue
        if g != g or age != age or g <= 0 or age < 50:
            continue
        rows.append((g, int(age), "men" if male == 1 else "women"))

def q(v, p):
    v = sorted(v); n = len(v)
    if n < 20: return None
    h = (n - 1) * p / 100.0
    lo = int(math.floor(h)); hi = min(lo + 1, n - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])

print("%-6s %-7s %5s %6s %8s %8s %9s" % ("sex","band","n_CH","P_cut_CH","P_cut_SH","diff","prev_CH"))
out = []
for sex in ("men","women"):
    for b in BANDS:
        ch = [g for g, a, s in rows if s == sex and b[0] <= a <= b[1]]
        if len(ch) < 30:
            continue
        cut = CUT[sex]
        pch = 100.0 * sum(1 for g in ch if g < cut) / len(ch)
        sh_vals = SHARE[(sex, b)]
        psh = pctl_pos(sh_vals, PCTL, cut)
        dif = (pch - psh) if psh is not None else None
        print("%-6s %-7s %5d %7.1f %8s %8s %8.1f" % (
            sex, "%d-%d" % (b[0], min(b[1],90)), len(ch), pch,
            ("<5" if psh is None else "%.1f" % psh),
            ("n/a" if dif is None else "%+.1f" % dif), pch))
        out.append((sex, b, len(ch), pch, psh, dif))

print("\n--- pooled 50+ ---")
for sex in ("men","women"):
    ch = [g for g, a, s in rows if s == sex]
    cut = CUT[sex]
    print(sex, "n=", len(ch), " prev=%.1f%%" % (100*sum(1 for g in ch if g < cut)/len(ch)),
          " P50=%.0f" % q(ch,50), " P25=%.0f" % q(ch,25), " P75=%.0f" % q(ch,75))

print("\n--- CHARLS percentile table for the report (men) ---")
for b in BANDS:
    ch = [g for g,a,s in rows if s=="men" and b[0]<=a<=b[1]]
    if len(ch) < 30: continue
    print("men %d-%d n=%d " % (b[0], min(b[1],90), len(ch)) + " ".join("%.0f" % q(ch,p) for p in PCTL))
print("--- women ---")
for b in BANDS:
    ch = [g for g,a,s in rows if s=="women" and b[0]<=a<=b[1]]
    if len(ch) < 30: continue
    print("women %d-%d n=%d " % (b[0], min(b[1],90), len(ch)) + " ".join("%.0f" % q(ch,p) for p in PCTL))
