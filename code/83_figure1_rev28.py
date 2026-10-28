# -*- coding: utf-8 -*-
"""83_figure1_rev28.py -- Figure 1 (cohort selection flow) as submitted.

Run from the repository root:  python3 code/83_figure1_rev28.py

This is the layout script behind the submitted Figure 1. Against the earlier
draft figure it changes typography only:

  * font DejaVu Sans -> Times New Roman, 6.8-8.2 pt -> 8.0-8.8 pt
  * a dashed leader is drawn from the main chain to the "not measured at both
    waves" box, which was unconnected in the draft
  * the canvas is cropped to the content, removing a 56.5 pt / 34.0 pt left-right
    asymmetry

Box geometry, connectors, colours and every number are unchanged: they are the
counts reported in Section 2.1 and in `manifest/outputs.csv` (Figure 1).

Requires matplotlib. Not part of `run_all.sh`: the pipeline's figure step is R,
and this script exists so the submitted layout can be regenerated exactly.
Writes figures/Figure_1.pdf and figures/Figure_1.png (1000 dpi).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
import matplotlib.font_manager as fm

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "pdf.fonttype": 42,
    "lines.antialiased": True,
})

INK = "#1B2A41"
MAIN_F, MAIN_E = "#DCE6F1", "#7F9DB9"
EXC_F,  EXC_E  = "#FBE3DE", "#D08479"
WHE_F,  WHE_E  = "#FFF3D6", "#D8B25E"
GREY_F, GREY_E = "#F2F2F2", "#9A9A9A"

# ---- 原始几何（单位 pt，取自提交版矢量 PDF 的 drawing 提取）----
CONTENT = (33.3, 9.8, 470.4, 530.1)      # x0, y0, x1, y1
MARGIN = 14.0

BOXES = [
    # x0,     y0,    x1,     y1,    fill,   edge,  dashed, fs,   text
    (80.8,   9.8,  299.4,  48.9, MAIN_F, MAIN_E, False, 8.8,
     "CHARLS 2011 national baseline\nn = 17,708"),
    (80.8,  87.9,  299.4, 127.0, MAIN_F, MAIN_E, False, 8.8,
     "Valid grip measurement in 2011 (n = 13,579)\nand in 2013 (n = 12,861)"),
    (80.8, 166.1,  299.4, 205.2, MAIN_F, MAIN_E, False, 8.8,
     "Measured at both waves\nn = 9,125"),
    (80.8, 244.3,  299.4, 283.4, MAIN_F, MAIN_E, False, 8.8,
     "ADL-free at the 2013 landmark\nn = 6,070"),
    (80.8, 322.5,  299.4, 361.5, MAIN_F, MAIN_E, False, 8.8,
     "Analytic sample\nn = 6,052"),
    (299.4, 200.3, 470.4, 249.2, EXC_F,  EXC_E,  False, 8.2,
     "Excluded\nADL disability, n = 438\nADL data missing, n = 2,617"),
    (299.4, 278.5, 470.4, 327.3, EXC_F,  EXC_E,  False, 8.2,
     "Excluded\nincomplete exposure or\nage/sex data, n = 18"),
    (54.6,  400.6, 154.4, 439.7, WHE_F,  WHE_E,  False, 8.8,
     "Died\nn = 286"),
    (161.6, 400.6, 285.1, 439.7, WHE_F,  WHE_E,  False, 8.8,
     "ADL assessed\nn = 4,453"),
    (332.6, 400.6, 446.7, 439.7, GREY_F, GREY_E, True,  8.8,
     "No 2018 outcome\nn = 1,313 (excluded)"),
    (33.3,  481.2, 308.9, 530.1, MAIN_F, MAIN_E, False, 8.8,
     "Analytic cohort with observed five-year outcome\nn = 4,739\n"
     "never low 3,833; remitted 325; incident 376; persistent 205"),
    (304.1, 51.3,  465.7, 114.8, GREY_F, GREY_E, True,  8.0,
     "Not measured at\nboth waves\n2011 only, n = 4,454\n2013 only, n = 3,736\nneither, n = 393"),
]

# 主链竖直箭头：(x, y_start, y_end)
STEMS = [(190.1, 48.9, 87.9), (190.1, 127.0, 166.1), (190.1, 205.2, 244.3),
         (190.1, 283.4, 322.5), (190.1, 361.5, 381.1)]
# 向右的排除箭头：(y, x_start, x_end)
BRANCH_R = [(224.7, 190.1, 299.4), (302.9, 190.1, 299.4)]
# 底部：分流与汇合
SPLIT_Y = 381.1
SPLIT_XS = [104.5, 223.3, 389.7]
MERGE_Y = 459.3
MERGE_XS = [104.5, 223.3]
FINAL_ARROW = (163.9, 459.3, 481.2)

LW = 0.9
HEAD_L, HEAD_W = 3.6, 1.8


def build():
    x0, y0, x1, y1 = CONTENT
    W = (x1 - x0) + 2 * MARGIN
    H = (y1 - y0) + 2 * MARGIN
    fig = plt.figure(figsize=(W / 72, H / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0 - MARGIN, x1 + MARGIN)
    ax.set_ylim(y1 + MARGIN, y0 - MARGIN)   # y 向下
    ax.axis("off")

    def seg(p0, p1, dashed=False):
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=INK, lw=LW,
                solid_capstyle="butt", zorder=3,
                linestyle=(0, (2.4, 1.6)) if dashed else "-")

    def head(tip, direction):
        dx, dy = direction
        L = (dx * dx + dy * dy) ** 0.5
        ux, uy = dx / L, dy / L
        bx, by = tip[0] - ux * HEAD_L, tip[1] - uy * HEAD_L
        px, py = -uy, ux
        ax.add_patch(Polygon([[tip[0], tip[1]],
                              [bx + px * HEAD_W, by + py * HEAD_W],
                              [bx - px * HEAD_W, by - py * HEAD_W]],
                             closed=True, fc=INK, ec=INK, lw=0.5, zorder=4))

    def arrow(p0, p1, dashed=False):
        seg(p0, p1)
        head(p1, (p1[0] - p0[0], p1[1] - p0[1]))

    # 框
    for bx0, by0, bx1, by1, fc, ec, dashed, fs, text in BOXES:
        ax.add_patch(FancyBboxPatch((bx0, by0), bx1 - bx0, by1 - by0,
                                    boxstyle="round,pad=0,rounding_size=5.5",
                                    lw=LW, fc=fc, ec=ec, zorder=2,
                                    linestyle=(0, (2.4, 1.6)) if dashed else "-"))
        ax.text((bx0 + bx1) / 2, (by0 + by1) / 2, text, ha="center", va="center",
                fontsize=fs, color=INK, linespacing=1.42, zorder=5)

    # 主链
    for x, ya, yb in STEMS:
        arrow((x, ya), (x, yb))
    # 排除分支
    for y, xa, xb in BRANCH_R:
        arrow((xa, y), (xb, y))
    # 底部三分流
    seg((190.1, 361.5), (190.1, SPLIT_Y))
    for xs in SPLIT_XS:
        seg((190.1, SPLIT_Y), (xs, SPLIT_Y))
        arrow((xs, SPLIT_Y), (xs, 400.6))
    # 汇合
    for xs in MERGE_XS:
        seg((xs, 439.7), (xs, MERGE_Y))
    seg((MERGE_XS[0], MERGE_Y), (MERGE_XS[1], MERGE_Y))
    arrow((FINAL_ARROW[0], MERGE_Y), (FINAL_ARROW[0], FINAL_ARROW[2]))

    # 新增：17,708 -> 9,125 这一段右侧的"未两次测量"虚线引线
    seg((190.1, 68.5), (295.0, 68.5), dashed=True)
    head((304.1, 68.5), (1, 0))
    return fig, ax


def fit_check(fig, ax):
    """框内文字是否溢出框"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bad = 0
    for (bx0, by0, bx1, by1, *_rest) in [(b[0], b[1], b[2], b[3]) for b in BOXES]:
        pass
    for b in BOXES:
        bx0, by0, bx1, by1, fc, ec, dashed, fs, text = b
        t = [c for c in ax.texts if c.get_text() == text][0]
        tb = t.get_window_extent(r)
        inv = ax.transData.inverted()
        (dx0, dy0) = inv.transform((tb.x0, tb.y0))
        (dx1, dy1) = inv.transform((tb.x1, tb.y1))
        over = max(bx0 - dx0, dx1 - bx1, by0 - min(dy0, dy1), max(dy0, dy1) - by1)
        if over > 0.5:
            bad += 1
            print("  OVERFLOW %5.2f pt | %s" % (over, text.split("\n")[0][:40]))
    print("fit check: %d boxes, %d overflow" % (len(BOXES), bad))
    return bad


if __name__ == "__main__":
    fig, ax = build()
    fit_check(fig, ax)
    fig.savefig(os.path.join(OUT, "Figure_1.pdf"), format="pdf")
    fig.savefig(os.path.join(OUT, "Figure_1.png"), dpi=1000, facecolor="white")
    print("written figures/Figure_1.pdf and figures/Figure_1.png")
