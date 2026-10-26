# -*- coding: utf-8 -*-
"""Build Supplementary Tables S1-S3 docx from the analysis CSVs."""
import csv, json

OUT = "out/Supplementary_Tables_S1-S3.docx"

def rd(p):
    with open(p, encoding="utf8") as f:
        return list(csv.DictReader(f))

rrr = rd("out/S2_weighted_MI_rrr.csv")      # model = "age+sex,19" / "fully adj,19"
absr = rd("out/S2_weighted_MI_absolute_risk.csv")
des = rd("out/S2_design_binomial.csv")
unw = rd("out/S2_multinom_weighted.csv")    # model = "weighted + cluster-robust (19)"
s3 = rd("out/S3_incremental.csv")

def f2(x): return f"{float(x):.2f}"
def ci(r): return f"{f2(r['RRR'])} ({f2(r['lo'])}-{f2(r['hi'])})"

def pick(rows, **kw):
    return [r for r in rows if all(str(r.get(k, "")).endswith(v) if k == "model_ends" else r.get(k) == v for k, v in kw.items())]

def rrr_table(cut):
    lines = ["| Outcome | Group | Age- and sex-adjusted | Fully adjusted |", "| --- | --- | --- | --- |"]
    for oc, lab in [("disabled", "ADL disability"), ("dead", "Death")]:
        for g in ["remitted", "incident", "persistent"]:
            a = [r for r in rrr if r["model"] == f"age+sex,{cut}" and r["outcome"] == oc and r["group"] == g][0]
            b = [r for r in rrr if r["model"] == f"fully adj,{cut}" and r["outcome"] == oc and r["group"] == g][0]
            lines.append(f"| {lab} | {g.capitalize()} | {ci(a)} | {ci(b)} |")
    return "\n".join(lines)

def abs_rows(cut, field):
    out = []
    for a in ["60", "65", "70", "75", "80"]:
        row = [r for r in absr if r["cutoff"] == cut and r["age"] == a]
        vals = [next(y for y in row if y["group"] == g) for g in ["never", "remitted", "incident", "persistent"]]
        out.append("| " + a + " | " + " | ".join(f"{float(x[field]):.1f}" for x in vals) + " |")
    return out

def abs_table():
    L = ["| Age | Never low | Remitted | Incident | Persistent |", "| --- | --- | --- | --- | --- |"]
    for title, cut, field in [("ADL disability, 2019 cut-offs", "19", "disabled"),
                              ("ADL disability, 2025 cut-offs", "25", "disabled"),
                              ("Death, 2019 cut-offs", "19", "dead"),
                              ("Death, 2025 cut-offs", "25", "dead")]:
        L.append(f"| **{title}** | | | | |")
        L += abs_rows(cut, field)
    return "\n".join(L)

def s2_table(cut):
    L = ["| Outcome | Group | Unweighted (cluster-robust) | Weighted (cluster-robust) | Design-based binomial |",
         "| --- | --- | --- | --- | --- |"]
    for oc, lab in [("disabled", "ADL disability"), ("dead", "Death")]:
        con = "disabled vs alive-free" if oc == "disabled" else "dead vs alive-free"
        for g in ["remitted", "incident", "persistent"]:
            u = [r for r in unw if r["model"].startswith("unweighted") and r["model"].endswith(f"({cut})") and r["outcome"] == oc and r["group"] == g][0]
            w = [r for r in unw if r["model"].startswith("weighted") and r["model"].endswith(f"({cut})") and r["outcome"] == oc and r["group"] == g][0]
            d = [r for r in des if r["cutoff"] == cut and r["contrast"] == con and r["group"] == g][0]
            L.append(f"| {lab} | {g.capitalize()} | {ci(u)} | {ci(w)} | {ci(d)} |")
    return "\n".join(L)

def s3_table():
    d = [r for r in s3 if r["outcome"].startswith("disability")][0]
    x = [r for r in s3 if r["outcome"].startswith("death")][0]
    L = ["| Metric | ADL disability (391 events) | Death (237 events) |", "| --- | --- | --- |"]
    def row(lbl, k, fmt=lambda v: f"{float(v):.3f}"):
        L.append(f"| {lbl} | {fmt(d[k])} | {fmt(x[k])} |")
    row("AUC, 2013 grip only (A)", "auc_A")
    row("AUC, A + 2011 low-grip flag (B1)", "auc_B1")
    row("AUC, A + four-group pattern (B2)", "auc_B2")
    row("DeLong P, A vs B1", "p_delong_B1")
    row("DeLong P, A vs B2", "p_delong_B2")
    row("Change in AIC, B1 vs A", "dAIC_B1", lambda v: f"{float(v):+.2f}")
    row("Change in AIC, B2 vs A", "dAIC_B2", lambda v: f"{float(v):+.2f}")
    row("Likelihood ratio chi-square, B1 vs A", "LR_B1", lambda v: f"{float(v):.2f}")
    row("P, B1 vs A", "p_LR_B1")
    row("Likelihood ratio chi-square, B2 vs A", "LR_B2", lambda v: f"{float(v):.2f}")
    row("P, B2 vs A", "p_LR_B2")
    row("Continuous NRI, B1 vs A", "NRI_B1", lambda v: f"{float(v):+.3f}")
    row("Continuous NRI, B2 vs A", "NRI_B2", lambda v: f"{float(v):+.3f}")
    row("Continuous NRI, B2 vs A, 95% CI lower", "NRI_B2_lo", lambda v: f"{float(v):+.3f}")
    row("Continuous NRI, B2 vs A, 95% CI upper", "NRI_B2_hi", lambda v: f"{float(v):+.3f}")
    row("IDI, B2 vs A", "IDI_B2", lambda v: f"{float(v):+.4f}")
    return "\n".join(L)

md = f"""# Supplementary Material

**Transient Versus Persistent Low Grip Strength and Five-Year Disability and Mortality in Middle-Aged and Older Chinese Adults**

## Supplementary Table S1. Re-analysis with the AWGS 2025 age-specific thresholds

**Panel A. Reclassification.** Low grip was defined as <28 kg in men and <18 kg in women at all ages (AWGS 2019), or as <34 kg in men and <20 kg in women at ages 50 to 64 years with <28 kg and <18 kg at ages 65 years and older (AWGS 2025). All rows are computed in the analytic sample (n = 6,052; recomputed from out/jamda_full.rds, R 4.5.2, 2026-10-23); among the 9,116 participants with classifiable grip measurements at both waves, the reversion proportion under the 2019 cut-offs was 59.8%.

| Characteristic | AWGS 2019 | AWGS 2025 |
| --- | --- | --- |
| Never low, n | 4,943 | 4,476 |
| Remitted, n | 385 | 575 |
| Incident, n | 462 | 619 |
| Persistent, n | 262 | 382 |
| Newly classified as low in 2013, n | - | 277 |
| Newly classified as low in 2011, n | - | 310 |
| Re-labelled as normal, n | 0 | 0 |
| Observed agreement, 2011 vs 2013 | 0.860 | 0.803 |
| Cohen's kappa | 0.304 | 0.273 |
| Low grip in 2011 that reverted by 2013, % | 59.5 | 60.1 |

**Panel B. Relative risk ratios (95% CI), weighted multiply imputed multinomial model, AWGS 2019 cut-offs.**

{rrr_table("19")}

*AWGS 2025 cut-offs (sensitivity analysis):*

{rrr_table("25")}

**Panel C. Standardized five-year risks (%), weighted multiply imputed model.**

{abs_table()}

Reference group for all relative risk ratios: never low. Models are weighted with the 2013 individual sampling weight, standard errors are clustered by community (408 clusters), and missing covariates are multiply imputed (20 datasets); covariate and outcome definitions follow the Methods of the main manuscript.

## Supplementary Table S2. Sampling weights, clustering and design-based inference

**Panel A. AWGS 2019 cut-offs.**

{s2_table("19")}

**Panel B. AWGS 2025 cut-offs.**

{s2_table("25")}

**Panel C. Choice of sampling weight (complete-case, fully adjusted multinomial models with cluster-robust standard errors, AWGS 2019 cut-offs, point estimates).**

| Weight | n | Disability: remitted / incident / persistent | Death: remitted / incident / persistent |
| --- | --- | --- | --- |
| 2011 individual weight | 4,459 | 1.29 / 1.56 / 1.95 | 0.88 / 1.30 / 1.88 |
| 2013 individual weight (used in the primary analysis) | 4,366 | 1.21 / 1.53 / 1.83 | 0.86 / 1.28 / 1.88 |
| 2013 longitudinal weight | 4,459 | 1.29 / 1.56 / 1.95 | 0.86 / 1.26 / 1.88 |

The 2018 weights cannot be used here: participants who died before 2018 have no 2018 weight, so weighting with them would remove all 286 deaths. Design-based models use community as the primary sampling unit, 52 strata defined by province and urban or rural residence, and the 2013 individual weight, with adjustment for single-PSU strata; the design degrees of freedom were 356. Panel C shows complete-case (non-imputed) estimates and does not reproduce the multiply imputed primary analysis; with the 2011 and 2013 longitudinal weights the persistent-group disability estimate was 1.95, above the 1.72 to 1.87 range of the multiply imputed specifications reported in the main text.

## Supplementary Table S4. Agreement between the 2011 and 2013 low-grip classifications (analytic sample, n = 6,052)

**Panel A. Two-by-two transition table.** Classifications use the AWGS 2019 cut-offs (<28 kg in men and <18 kg in women) applied at all ages. The four cells are the four exposure groups of the main analysis: persistent (low at both waves), remitted (low in 2011, normal in 2013), incident (normal in 2011, low in 2013), and never low (normal at both waves).

| | Low grip in 2013 | Normal grip in 2013 | Total |
| --- | --- | --- | --- |
| **Low grip in 2011** | 262 | 385 | 647 |
| **Normal grip in 2011** | 462 | 4,943 | 5,405 |
| **Total** | 724 | 5,328 | 6,052 |

Cohen's kappa = 0.30; observed agreement 86.0% (5,205/6,052), expected agreement 79.9%. Denominator is the analytic sample (n = 6,052); transition proportions reported in the main text are based on all 9,125 participants measured at both waves (59.8% of the 859 with low grip in 2011 were normal in 2013) and on the 5,405 analytic-sample participants with normal grip in 2011 (8.6% became low), respectively.

## Supplementary Table S3. Incremental predictive value of the 2011 reading

Model A = covariates + 2013 grip strength (per 5 kg). Model B1 = A + 2011 low-grip flag. Model B2 = A + four-group pattern. Binary outcomes; n = 4,222 for disability and 4,068 for death.

{s3_table()}

With the AWGS 2025 cut-offs for the four-group pattern the results were unchanged: AUC 0.695 to 0.696 for disability (DeLong P = 0.155; change in AIC +4.81; likelihood ratio P = 0.755) and 0.785 to 0.785 for death (DeLong P = 0.68; change in AIC +5.75; likelihood ratio P = 0.97). AUC, area under the receiver operating characteristic curve; NRI, net reclassification improvement; IDI, integrated discrimination improvement.
"""
open("_supp.md", "w", encoding="utf8").write(md)
print("markdown written,", len(md), "chars")

# write the assembled tables to a plain .docx (no external editor needed)
try:
    import docx
    d = docx.Document()
    for block in md.split("\n\n"):
        block = block.strip()
        if block:
            d.add_paragraph(block)
    d.save(OUT)
    print("docx written:", OUT)
except Exception as e:
    print("docx step skipped:", e)
