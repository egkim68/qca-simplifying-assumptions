# figures.R -- Figures 1-4 of the paper from results/. Usage: Rscript figures.R
d <- "results/"; o <- "figures/"; dir.create(o, showWarnings = FALSE)
cl <- read.csv(paste0(d, "figure2.csv"))
pv <- read.csv(paste0(d, "figure3.csv"))
sg <- read.csv(paste0(d, "figure4a.csv"))
nz <- c(0, 0.05, 0.1, 0.2); nzl <- c("0.00", "0.05", "0.10", "0.20")
GREEN <- "#1E5B47"; BLUE <- "#2A5F8F"; LGREEN <- "#A9CDBF"; LBLUE <- "#A8C3DD"; MID <- "#5E9C86"
W <- 6.5
dev0 <- function(f, h) png(paste0(o, f), width = W, height = h, units = "in", res = 600, pointsize = 10)
par0 <- function() par(family = "sans", mar = c(3.6, 4.2, 2.2, 0.4), mgp = c(2.4, 0.6, 0), las = 1, bty = "l",
                       cex.axis = 0.95, cex.lab = 1.05, cex.main = 1.05, font.main = 1, tcl = -0.3)

# Figure 1: workflow
dev0("fig1.png", 1.8)
par(mar = c(0, 0, 0, 0), family = "sans"); plot.new(); plot.window(c(0, 10), c(0, 2))
lab <- c("True Boolean\nfunction", "Sample cases\n(with or without\nnoise)", "Truth table\n(main cut-off\n0.80)", "Solution\nand checks", "Score against\nthe truth")
x <- seq(1, 9, length.out = 5)
for (i in 1:5) {
  rect(x[i] - 0.8, 0.3, x[i] + 0.8, 1.7, col = if (i %% 2) LGREEN else LBLUE, border = if (i %% 2) GREEN else BLUE, lwd = 2)
  text(x[i], 1, lab[i], cex = 0.92, col = "black")
}
for (i in 1:4) arrows(x[i] + 0.82, 1, x[i + 1] - 0.82, 1, length = 0.08, lwd = 2, col = GREEN)
dev.off()

# Figure 2: cell-level gaps
dev0("fig2.png", 3.4); par(mfrow = c(1, 2)); par0()
boxplot(h1 ~ noise, cl, col = LGREEN, border = GREEN, outline = FALSE, lwd = 1.3, names = nzl, xlab = "Outcome noise rate", ylab = "Precision gap", main = "(a) Non-CSA vs CSA")
abline(h = 0.15, lty = 2, lwd = 1.5, col = BLUE); abline(h = 0, col = "grey60")
boxplot(h2 ~ noise, cl, col = LBLUE, border = BLUE, outline = FALSE, lwd = 1.3, names = nzl, xlab = "Outcome noise rate", ylab = "Precision gap", main = "(b) Tenable vs untenable")
abline(h = 0.10, lty = 2, lwd = 1.5, col = GREEN); abline(h = 0, col = "grey60")
dev.off()

# Figure 3: parsimonious vs cleaned
dev0("fig3.png", 3.4); par(mfrow = c(1, 2)); par0()
for (v in c("sufficient", "bt_correct")) {
  a <- pv[pv$solution == "pars", ]; a <- a[order(a$noise), ]; b <- pv[pv$solution == "clean", ]; b <- b[order(b$noise), ]
  barplot(rbind(a[[v]], b[[v]]), beside = TRUE, names.arg = nzl, col = c(GREEN, LGREEN), border = GREEN, lwd = 1.3, ylim = c(0, 1),
          xlab = "Outcome noise rate", ylab = "Share of terms", main = ifelse(v == "sufficient", "(a) Truly sufficient", "(b) Meets B&T criterion"))
  if (v == "sufficient") legend("topright", c("Parsimonious", "Cleaned"), fill = c(GREEN, LGREEN), border = GREEN, bty = "n", cex = 1.0)
}
dev.off()

# Figure 4: support and AUC
dev0("fig4.png", 3.9); par(mfrow = c(1, 2)); par0(); par(cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.15)
gs <- c("1-2", "3-4", "5+"); cols <- c(BLUE, MID, GREEN); pc <- c(21, 22, 24)
plot(NA, xlim = range(nz), ylim = c(0, 1), xaxt = "n", xlab = "Outcome noise rate", ylab = "Share truly sufficient", main = "(a) By supporting rows"); axis(1, nz, nzl)
for (i in 1:3) { s <- sg[sg$sgrp == gs[i], ]; s <- s[order(s$noise), ]
  lines(s$noise, s$sufficient, type = "b", pch = pc[i], bg = cols[i], col = cols[i], lwd = 2, cex = 1.2) }
legend("bottomleft", paste(gs, "rows"), pch = pc, pt.bg = cols, col = cols, lwd = 2, bty = "n", cex = 1.15)
rv <- read.csv(paste0(d, "figure4b.csv")); rv <- rv[order(rv$noise), ]
plot(rv$noise, rv$full, type = "b", pch = 21, bg = GREEN, col = GREEN, lwd = 2, cex = 1.1, ylim = c(0.4, 0.9), xaxt = "n",
     xlab = "Outcome noise rate", ylab = "AUC (test half)", main = "(b) Predicting sufficiency"); axis(1, nz, nzl)
lines(rv$noise, rv$rem_share, type = "b", pch = 24, bg = MID, col = MID, lwd = 2, cex = 1.1)
lines(rv$noise, rv$support, type = "b", pch = 23, bg = LBLUE, col = BLUE, lwd = 2, cex = 1.1, lty = 3)
lines(rv$noise, rv$raw_cov, type = "b", pch = 22, bg = BLUE, col = BLUE, lwd = 2, cex = 1.1, lty = 2)
legend("bottomleft", c("Support model", "Remainder share only", "Support count only", "Raw coverage"), pch = c(21, 24, 23, 22),
       pt.bg = c(GREEN, MID, LBLUE, BLUE), col = c(GREEN, MID, BLUE, BLUE), lwd = 2, lty = c(1, 1, 3, 2), bty = "n", cex = 1.0)
dev.off()
cat("done\n")
