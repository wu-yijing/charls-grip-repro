# 07b_figure2.R -- Figure 2: (A) observed outcome distribution, (B) observed vs expected transitions
dir.create("figures", showWarnings = FALSE)

BLUE   <- "#3B6EA5"
ORANGE <- "#D08B4F"
GREY   <- "#9A9A9A"
LIGHT  <- "#C9D6E8"

png("figures/Figure_2.png", width = 3640, height = 3820,
    res = 400, type = "cairo", family = "serif")

layout(matrix(1:2, nrow = 2), heights = c(1.0, 1.05))

## ---------------- Panel A ----------------
par(mar = c(5.4, 5.6, 4.6, 1.4), mgp = c(3.4, 0.8, 0), cex.axis = 1.05, cex.lab = 1.18, xpd = NA)

m <- matrix(c(87.8, 79.1, 73.4, 62.0,
              7.7, 11.7, 15.2, 17.1,
              4.4,  9.2, 11.4, 21.0), nrow = 3, byrow = TRUE)
nm <- c("Never low\n(n = 3,833)", "Remitted\n(n = 325)",
        "Incident\n(n = 376)", "Persistent\n(n = 205)")

bp <- barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, ylim = c(0, 112),
              ylab = "Participants, %", names.arg = nm, las = 1,
              space = 0.5, axes = FALSE, cex.names = 1.0)
axis(2, at = seq(0, 100, 20), las = 1, tick = FALSE, col.axis = "grey25")
abline(h = seq(0, 100, 20), col = "grey88", lwd = 0.8)
barplot(m, col = c(BLUE, ORANGE, GREY), border = NA, add = TRUE, axes = FALSE, space = 0.5)

## segment labels omitted: the numbers are given in Table 2
legend(x = mean(bp), y = 110, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n", cex = 1.05,
       fill = c(BLUE, ORANGE, GREY), border = NA,
       legend = c("Alive, ADL-free", "ADL disability", "Death"))
mtext("A", side = 3, line = 2.6, adj = -0.13, font = 2, cex = 1.5)

## ---------------- Panel B ----------------
par(mar = c(8.2, 5.6, 4.0, 1.4))

lab <- c("Reversion\n(low in 2011, normal in 2013)", "Incidence\n(normal in 2011, low in 2013)")
obs <- c(59.5, 8.5);  exp <- c(55.6, 10.8);  lo <- c(52.5, 10.1);  hi <- c(58.9, 11.4)

b <- barplot(rbind(obs, exp), beside = TRUE, border = NA,
             col = c(BLUE, LIGHT), ylim = c(0, 88), ylab = "Participants, %",
             names.arg = c("", ""), space = c(0.3, 1.2), axes = FALSE)
axis(2, at = seq(0, 80, 10), las = 1, tick = FALSE, col.axis = "grey25")
abline(h = seq(0, 80, 10), col = "grey88", lwd = 0.8)
barplot(rbind(obs, exp), beside = TRUE, border = NA, col = c(BLUE, LIGHT),
        add = TRUE, axes = FALSE, space = c(0.3, 1.2))

arrows(b[2, ], lo, b[2, ], hi, angle = 90, code = 3, length = 0.045, lwd = 1.6, col = "grey25")
text(b[1, ], obs + 2.6, sprintf("%.1f", obs), cex = 1.0, col = "grey15")
text(b[2, ], hi + 2.6, sprintf("%.1f", exp), cex = 1.0, col = "grey15")

text(mean(b[, 1]), -8.0, "Reversion", xpd = NA, cex = 1.0)
text(mean(b[, 1]), -13.0, "(low in 2011, normal in 2013)", xpd = NA, cex = 0.90)
text(mean(b[, 2]), -8.0, "Incidence", xpd = NA, cex = 1.0)
text(mean(b[, 2]), -13.0, "(normal in 2011, low in 2013)", xpd = NA, cex = 0.90)

legend(x = max(b) + 0.15, y = 84, xjust = 1, yjust = 0.5, bty = "n", cex = 1.05,
       fill = c(BLUE, LIGHT), border = NA,
       legend = c("Observed", "Expected if muscle function did not change"))
mtext("B", side = 3, line = 2.6, adj = -0.13, font = 2, cex = 1.5)

invisible(dev.off())
cat("written\n")
