# 84_figure2_rev28.R -- Figure 2 as submitted (typography only; values unchanged).
#
# Run from the repository root:  Rscript code/84_figure2_rev28.R
#
# Against 81_figure2_rev23.R this changes layout only:
#   1) font sizes normalised to ~9 pt (they were 9.8-12.6 pt at 166 mm width)
#   2) panel A: ylim 96 -> 112 and the legend moved above the bars, so it no
#      longer sits on the 80% grid line
#   3) panel B: horizontal grid lines dropped; at 20% spacing they ran through the
#      data labels in 12 places (the submitted figure has none)
#   4) panel B: segments below 6% are labelled outside the bar, level with the
#      segment, instead of inside a 4.4%-tall block
#   5) top margin widened so the panel letters A/B do not touch the edge
# Values, colours, panel order and legend wording are unchanged.
#
# Requires R with cairo (cairo_pdf, png type="cairo"). Not part of run_all.sh:
# 07b_figure2.R and 81_figure2_rev23.R already produce Figure 2 in the pipeline;
# this script produces the submitted layout.
# Writes figures/Figure_2.pdf and figures/Figure_2.png (1000 dpi).
# The submitted file is this canvas trimmed to the ink box (470.0 x 467.2 pt =
# 165.8 x 164.8 mm, i.e. under the 166 mm text width); trimming removes outer
# white space only. All 52 text strings are identical to the submitted PDF.

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BLUE   <- "#3B6EA5"
ORANGE <- "#D08B4F"
GREY   <- "#9A9A9A"
LIGHT  <- "#C9D6E8"

W_IN <- 6703 / 1000
H_IN <- 6925 / 1000

draw <- function() {
  par(mfrow = c(2, 1), omi = c(0.15, 0, 0, 0), ps = 12)

  ## ---------------- Panel A ----------------
  par(mar = c(4.4, 5.0, 3.2, 1.2), mgp = c(2.9, 0.65, 0), family = "serif")

  obs <- c(59.5, 8.5); olo <- c(55.7, 7.8); ohi <- c(63.3, 9.2)
  exp <- c(55.6, 10.8); elo <- c(52.5, 10.1); ehi <- c(58.9, 11.4)

  b <- barplot(rbind(obs, exp), beside = TRUE, border = NA, col = c(BLUE, LIGHT),
               ylim = c(0, 112), ylab = "Proportion, %", axes = FALSE,
               space = c(0.12, 1.15), width = 0.42, names.arg = c("", ""),
               cex.axis = 0.75, cex.lab = 0.79)
  axis(2, at = seq(0, 80, 10), las = 1, tick = TRUE, tcl = -0.12, lwd.ticks = 0.7,
       col.axis = "grey25", cex.axis = 0.75)
  # 去掉横向网格线：柱状图的刻度已足够读数，网格线是"线条穿字"的根源
  barplot(rbind(obs, exp), beside = TRUE, border = NA, col = c(BLUE, LIGHT),
          add = TRUE, axes = FALSE, space = c(0.12, 1.15), width = 0.42)

  arrows(b[1, ], olo, b[1, ], ohi, angle = 90, code = 3, length = 0.05, lwd = 1.7, col = "grey20")
  arrows(b[2, ], elo, b[2, ], ehi, angle = 90, code = 3, length = 0.05, lwd = 1.7, col = "grey35")
  text(b[1, ], ohi, sprintf("%.1f", obs), pos = 3, offset = 0.55, cex = 0.71, col = "grey10")
  text(b[2, ], ehi, sprintf("%.1f", exp), pos = 3, offset = 0.55, cex = 0.71, col = "grey10")

  mtext("Reversion", side = 1, line = 1.9, at = mean(b[, 1]), cex = 0.75)
  mtext("(low in 2011, normal in 2013)", side = 1, line = 3.2, at = mean(b[, 1]), cex = 0.68)
  mtext("Incidence", side = 1, line = 1.9, at = mean(b[, 2]), cex = 0.75)
  mtext("(normal in 2011, low in 2013)", side = 1, line = 3.2, at = mean(b[, 2]), cex = 0.68)

  # 图例放在绘图区顶部（ylim=112 -> 图例下沿约在 100，高于 80 网格线）
  legend("topright", bty = "n", cex = 0.71, fill = c(BLUE, LIGHT), border = NA,
         legend = c("Observed", "Expected, no individual change (2013 marginal)"))
  mtext("A", side = 3, line = 1.4, adj = -0.09, font = 2, cex = 1.0)

  ## ---------------- Panel B ----------------
  par(mar = c(5.2, 5.0, 3.4, 1.2))
  m <- matrix(c(87.8, 79.1, 73.4, 62.0,
                7.7, 11.7, 15.2, 17.1,
                4.4,  9.2, 11.4, 21.0), nrow = 3, byrow = TRUE)
  nm <- c("Never low\n(n = 3,833)", "Remitted\n(n = 325)",
          "Incident\n(n = 376)", "Persistent\n(n = 205)")
  SP <- 0.55
  bp <- barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, ylim = c(0, 118),
                ylab = "Proportion, %", names.arg = rep("", 4), las = 1, space = SP,
                axes = FALSE, cex.axis = 0.75, cex.lab = 0.79)
  mtext(sub("\n.*", "", nm), side = 1, at = bp, line = 0.65, cex = 0.75)
  mtext(sub(".*\n", "", nm), side = 1, at = bp, line = 1.8, cex = 0.68)
  axis(2, at = seq(0, 100, 20), las = 1, tick = TRUE, tcl = -0.12, lwd.ticks = 0.7,
       col.axis = "grey25", cex.axis = 0.75)
  # 同面板 A：不画横向网格线
  barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, add = TRUE, axes = FALSE, space = SP)

  halfw <- 0.5   # barplot 中单柱宽度为 1，半宽即 0.5；space 是间隙
  ymid <- apply(m, 2, cumsum) - m / 2
  for (j in 1:4) for (i in 1:3) {
    lab <- sprintf("%.1f", m[i, j])
    if (m[i, j] < 6) {   # 只有真正放不下的薄段（4.4%）才外置
      # 段太薄，标注放到柱右侧，垂直对齐该段中线
      text(bp[j] + halfw + 0.08, ymid[i, j], lab, adj = c(0, 0.5),
           cex = 0.71, col = "grey20")
    } else {
      text(bp[j], ymid[i, j], lab, cex = 0.71,
           col = if (i == 1) "white" else "grey10")
    }
  }

  legend(x = mean(bp), y = 112, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
         cex = 0.71, fill = c(BLUE, ORANGE, GREY), border = NA,
         legend = c("Alive, ADL-free", "ADL disability", "Death"))
  mtext("B", side = 3, line = 1.4, adj = -0.09, font = 2, cex = 1.0)
}

png(file.path(OUT, "Figure_2.png"), width = 6703, height = 6925, res = 1000,
    type = "cairo", family = "serif", bg = "white")
draw(); invisible(dev.off())

cairo_pdf(file.path(OUT, "Figure_2.pdf"), width = W_IN, height = H_IN, family = "serif")
draw(); invisible(dev.off())

cat("written figures/Figure_2.pdf and figures/Figure_2.png\n")
