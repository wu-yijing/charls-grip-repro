# -*- coding: utf-8 -*-
"""Figure 2 as submitted: observed reversion and incidence against the closed-form
no-change expectations. Every number is read from the aggregate inputs in
manifest/fig_inputs/ -- nothing is hard-coded and no respondent-level data are used.

  manifest/fig_inputs/JFA_boot_persist.csv          observed rates, cluster bootstrap (Table 1)
  manifest/fig_inputs/JFA_benchmark_reconcile.csv   closed-form expectations (Table S32, Panel B)

Output: figures/Figure_2.pdf and figures/Figure_2.png (RGB, 300 dpi).
Run from the repository root:  python code/85_figure2_submitted.py
"""
import os, csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IN = os.path.join(ROOT, "manifest", "fig_inputs")
FIG = os.path.join(ROOT, "figures")

boot = list(csv.DictReader(open(os.path.join(IN, "JFA_boot_persist.csv"))))
recon = list(csv.DictReader(open(os.path.join(IN, "JFA_benchmark_reconcile.csv"))))

def b(setname, metric, col):
    return float(next(r for r in boot if r["set"] == setname and r["metric"] == metric)[col])

obs = [100 * b("all_2019", "revert", "est"), 100 * b("all_2019", "incid", "est")]
obs_lo = [100 * b("all_2019", "revert", "lo"), 100 * b("all_2019", "incid", "lo")]
obs_hi = [100 * b("all_2019", "revert", "hi"), 100 * b("all_2019", "incid", "hi")]

order = ["Sex-specific reliability (men 0.66, women 0.51)",
         "Single reliability 0.674 (submitted specification)",
         "Single pooled reliability 0.723 (consistent with SD 9.71)"]
labels = ["No change: sex-specific reliability",
          "No change: reliability 0.674",
          "No change: pooled reliability 0.723"]
exp = {}
for key, lab in zip(order, labels):
    r = next(x for x in recon if x["spec"] == key)
    exp[lab] = [float(r["reversion"]), float(r["incidence"])]

bars = {"Observed": (obs, [obs_lo, obs_hi])}
for lab in labels:
    bars[lab] = (exp[lab], None)

x = np.arange(2); w = 0.19
fig, ax = plt.subplots(figsize=(9.2, 4.6))
shades = ["#3b3b3b", "#8c8c8c", "#b8b8b8", "#dcdcdc"]
for i, (name, (vals, ci)) in enumerate(bars.items()):
    off = (i - 1.5) * w
    err = None
    if ci:
        err = np.array([[vals[j] - ci[0][j] for j in range(2)],
                        [ci[1][j] - vals[j] for j in range(2)]])
    ax.bar(x + off, vals, w, label=name, color=shades[i], edgecolor="black", linewidth=0.6,
           yerr=err, capsize=3, error_kw={"linewidth": 0.8, "ecolor": "black"})
    if name != "Observed":
        for j, v in enumerate(vals):
            ax.text(x[j] + off, v + 0.9, "%.1f" % v, ha="center", va="bottom", fontsize=7)
ax.set_ylabel("Proportion, %")
ax.set_xticks(x); ax.set_xticklabels(["Reversion", "Incidence"])
ax.set_ylim(0, 78); ax.set_xlim(-0.55, 1.55)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="upper right")
for xi, v, hi in zip(x, obs, obs_hi):
    ax.text(xi - 1.5 * w, hi + 1.6, "%.1f" % v, ha="center", fontsize=8)
fig.tight_layout()

pdf = os.path.join(FIG, "Figure_2.pdf")
png = os.path.join(FIG, "Figure_2.png")
fig.savefig(pdf, format="pdf", dpi=300)
fig.savefig(png, format="png", dpi=300)
plt.close(fig)
Image.open(png).convert("RGB").save(png, dpi=(300, 300))
print("observed %.1f (%.1f-%.1f) / %.1f (%.1f-%.1f)" %
      (obs[0], obs_lo[0], obs_hi[0], obs[1], obs_lo[1], obs_hi[1]))
for lab in labels:
    print("%s: %.1f / %.1f" % (lab, exp[lab][0], exp[lab][1]))
print("mode", Image.open(png).mode, Image.open(png).size, Image.open(png).info.get("dpi"))
