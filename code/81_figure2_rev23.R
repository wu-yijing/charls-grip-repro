# Ported for the reproducibility package: run from the repository root.
# Draws Figure 2 (600 dpi PNG + vector PDF) from the values reported in Sections 3.3
# and 3.4 and Table 2. Writes to figures/.
# 81_figure2_rev23.R -- Figure 2, revised: the expected values are labelled for the
# calibration actually plotted (matched to the 2013 marginal, i.e. the observed mean
# decline is imposed; no *individual* change), and the whiskers are defined.

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE)

BLUE   <- "#3B6EA5"
ORANGE <- "#D08B4F"
GREY   <- "#9A9A9A"
LIGHT  <- "#C9D6E8"

W <- 4022; H <- 4155; RES <- 600

draw <- function() {
  par(mfrow = c(2, 1), omi = c(0.2, 0, 0, 0))

  ## ---------------- Panel A: observed vs expected ----------------
  par(mar = c(4.6, 5.4, 3.0, 1.4), mgp = c(3.1, 0.75, 0), cex.axis = 1.05,
      cex.lab = 1.16, family = "serif")

  obs <- c(59.5, 8.5); olo <- c(55.7, 7.8); ohi <- c(63.3, 9.2)
  exp <- c(55.6, 10.8); elo <- c(52.5, 10.1); ehi <- c(58.9, 11.4)

  b <- barplot(rbind(obs, exp), beside = TRUE, border = NA, col = c(BLUE, LIGHT),
               ylim = c(0, 96), ylab = "Proportion, %", axes = FALSE,
               space = c(0.12, 1.15), width = 0.42, names.arg = c("", ""))
  axis(2, at = seq(0, 80, 10), las = 1, tick = FALSE, col.axis = "grey25")
  abline(h = seq(0, 80, 10), col = "grey88", lwd = 0.9)
  barplot(rbind(obs, exp), beside = TRUE, border = NA, col = c(BLUE, LIGHT),
          add = TRUE, axes = FALSE, space = c(0.12, 1.15), width = 0.42)

  arrows(b[1, ], olo, b[1, ], ohi, angle = 90, code = 3, length = 0.05, lwd = 1.7, col = "grey20")
  arrows(b[2, ], elo, b[2, ], ehi, angle = 90, code = 3, length = 0.05, lwd = 1.7, col = "grey35")
  text(b[1, ], ohi, sprintf("%.1f", obs), pos = 3, offset = 0.55, cex = 1.05, col = "grey10")
  text(b[2, ], ehi, sprintf("%.1f", exp), pos = 3, offset = 0.55, cex = 1.05, col = "grey10")

  mtext("Reversion", side = 1, line = 2.0, at = mean(b[, 1]), cex = 1.05)
  mtext("(low in 2011, normal in 2013)", side = 1, line = 3.4, at = mean(b[, 1]), cex = 0.90)
  mtext("Incidence", side = 1, line = 2.0, at = mean(b[, 2]), cex = 1.05)
  mtext("(normal in 2011, low in 2013)", side = 1, line = 3.4, at = mean(b[, 2]), cex = 0.90)

  legend("topright", bty = "n", cex = 0.95, fill = c(BLUE, LIGHT), border = NA,
         legend = c("Observed", "Expected, no individual change (2013 marginal)"))
  mtext("A", side = 3, line = 1.5, adj = -0.10, font = 2, cex = 1.45)

  ## ---------------- Panel B: outcome distribution ----------------
  par(mar = c(5.8, 5.4, 3.0, 1.4))
  m <- matrix(c(87.8, 79.1, 73.4, 62.0,
                7.7, 11.7, 15.2, 17.1,
                4.4,  9.2, 11.4, 21.0), nrow = 3, byrow = TRUE)
  nm <- c("Never low\n(n = 3,833)", "Remitted\n(n = 325)",
          "Incident\n(n = 376)", "Persistent\n(n = 205)")

  bp <- barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, ylim = c(0, 128),
                ylab = "Proportion, %", names.arg = rep("", 4), las = 1, space = 0.55,
                axes = FALSE)
  mtext(sub("\n.*", "", nm), side = 1, at = bp, line = 0.7, cex = 1.0)
  mtext(sub(".*\n", "", nm), side = 1, at = bp, line = 1.9, cex = 0.92)
  axis(2, at = seq(0, 100, 20), las = 1, tick = FALSE, col.axis = "grey25")
  abline(h = seq(0, 100, 20), col = "grey88", lwd = 0.9)
  barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, add = TRUE, axes = FALSE, space = 0.55)

  ymid <- apply(m, 2, cumsum) - m / 2
  for (j in 1:4) for (i in 1:3) {
    yy <- ymid[i, j]; cx <- 0.95
    if (m[i, j] < 8) { cx <- 0.82 }
    text(bp[j], yy, sprintf("%.1f", m[i, j]), cex = cx,
         col = if (i == 1) "white" else "grey10")
  }

  legend(x = mean(bp), y = 124, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
         cex = 1.02, fill = c(BLUE, ORANGE, GREY), border = NA,
         legend = c("Alive, ADL-free", "ADL disability", "Death"))
  mtext("B", side = 3, line = 1.5, adj = -0.10, font = 2, cex = 1.45)
}

png(file.path(OUT, "FIG2_rev23.png"), width = W, height = H, res = RES,
    type = "cairo", family = "serif", bg = "white")
draw(); invisible(dev.off())

cairo_pdf(file.path(OUT, "FIG2_rev23.pdf"), width = W / RES, height = H / RES, family = "serif")
draw(); invisible(dev.off())

cat("written\n")
