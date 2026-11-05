# -*- coding: utf-8 -*-
"""Figure 3 as submitted.

Panel A (decision arithmetic) restates the published counts -- 107 low at the first reading of
1,000 screened, 43 low again, 64 not -- so its numbers are fixed by the manuscript and are not
read from a file. Panel B reads the distance inputs in manifest/fig_inputs/:

  manifest/fig_inputs/JFA_distance_bins.csv    observed proportion below the cut-off in 2013,
                                               by distance of the 2011 reading (Table S36)
  manifest/fig_inputs/JFA_distance_curve.csv   the no-change model curve at the same distances

Output: figures/Figure_3.pdf and figures/Figure_3.png (RGB, 300 dpi).
Run from the repository root:  python code/86_figure3_submitted.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, os
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IN = os.path.join(ROOT, "manifest", "fig_inputs")
FIG = os.path.join(ROOT, "figures")

bins = list(csv.DictReader(open(os.path.join(IN, "JFA_distance_bins.csv"))))
curve = list(csv.DictReader(open(os.path.join(IN, "JFA_distance_curve.csv"))))
xs = [float(r["mid"]) for r in bins]; ps = [float(r["p"]) * 100 for r in bins]
los = [float(r["lo"]) * 100 for r in bins]; his = [float(r["hi"]) * 100 for r in bins]
cx = [float(r["d"]) for r in curve]; cy = [float(r["p"]) * 100 for r in curve]

fig = plt.figure(figsize=(7.1, 3.75), dpi=300)
axA = fig.add_axes([0.01, 0.02, 0.29, 0.92]); axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
NAVY = "#2F4A6B"; LIGHT = "#E8EEF6"; EDGE = "#2E5C8A"; GREY = "#5F5F5F"
def box(x, y, w, h, t, fc=LIGHT, ec=EDGE, tc="#1A1A1A", dashed=False, fs=6.2):
    axA.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.02",
        linewidth=0.8, edgecolor=ec, facecolor=fc, linestyle=(0, (3, 2)) if dashed else "solid"))
    axA.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs, color=tc, linespacing=1.5)
def arr(x1, y1, x2, y2):
    axA.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=5,
        linewidth=0.7, color="#333333", shrinkA=0, shrinkB=0))
box(0.12, 0.86, 0.76, 0.085, "1,000 screened")
arr(0.50, 0.86, 0.50, 0.815)
box(0.03, 0.68, 0.94, 0.12, "107 low at the first\nreading (10.7%)")
arr(0.50, 0.68, 0.50, 0.638)
axA.plot([0.25, 0.75], [0.638, 0.638], color="#333333", lw=0.7)
arr(0.25, 0.638, 0.25, 0.598); arr(0.75, 0.638, 0.75, 0.598)
box(0.02, 0.44, 0.46, 0.155, "43 low again\n(40.5%)", fc=NAVY, ec=NAVY, tc="white")
box(0.52, 0.44, 0.46, 0.155, "64 not low\n(59.5%)", fc="white", ec=GREY, dashed=True)
axA.text(0.75, 0.405, "label not\nrepeated", ha="center", va="top", fontsize=5.8, color=GREY, style="italic", linespacing=1.4)
axA.plot([0.02, 0.98], [0.33, 0.33], color="#C4C4C4", lw=0.6)
axA.text(0.02, 0.30, "Without a repeat all 107 are\nlabeled at the first reading; 43\nwould have been low two years later.",
         ha="left", va="top", fontsize=5.8, color="#1A1A1A", linespacing=1.6)

axB = fig.add_axes([0.395, 0.17, 0.575, 0.70])
axB.errorbar(xs, ps, yerr=[[p - l for p, l in zip(ps, los)], [h - p for p, h in zip(ps, his)]],
             fmt="o", ms=4, mfc="white", mec=NAVY, mew=1.0, ecolor=NAVY, elinewidth=0.8, capsize=2,
             label="Observed (95% CI)", zorder=3)
axB.plot(cx, cy, color="#B4562A", lw=1.6, label="No-change model", zorder=2)
axB.axvline(0, color="#888888", ls=(0, (3, 2)), lw=0.8)
axB.text(0.25, 78, "cut-off", fontsize=6, color="#666666")
axB.set_xlim(-12.5, 8.5); axB.set_ylim(0, 85)
axB.set_xlabel("2011 reading minus the sex-specific cut-off, kg", fontsize=7)
axB.set_ylabel("Below the cut-off in 2013, %", fontsize=7)
axB.tick_params(labelsize=6.5, width=0.7, length=2.5)
for s in ("top", "right"): axB.spines[s].set_visible(False)
for s in ("left", "bottom"): axB.spines[s].set_linewidth(0.7)
axB.legend(frameon=False, fontsize=6.3, loc="upper right", handlelength=1.6, borderpad=0.2)
axB.annotate("63% (model)", xy=(-10, 63), xytext=(-9.4, 75), ha="center", fontsize=6.2, color="#B4562A",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="#B4562A"))
axB.annotate("27% (model)", xy=(0, 26.6), xytext=(2.0, 42), ha="center", fontsize=6.2, color="#B4562A",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="#B4562A"))

fig.text(0.155, 0.985, "A. Decision arithmetic", ha="center", va="top", fontsize=7, fontweight="bold")
fig.text(0.683, 0.985, "B. Probability that a low reading persists, by distance to the cut-off", ha="center", va="top", fontsize=7, fontweight="bold")

out = os.path.join(FIG, "Figure_3")
fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.02)
fig.savefig(out + ".png", bbox_inches="tight", pad_inches=0.02, facecolor="white")
from PIL import Image
Image.open(out + ".png").convert("RGB").save(out + ".png", dpi=(300, 300))
print("saved", out + ".pdf", out + ".png")
