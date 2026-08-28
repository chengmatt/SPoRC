# Purpose: Render the schematic in the reference point uncertainty section of
#          vignettes/i_reference_points.Rmd, which shows where the sensitivity
#          d_j comes from. Changing a parameter deforms the criterion and leaves
#          a non-zero slope at the old solution; the new solution is however far
#          you have to walk to shed that slope. Both panels get the same
#          parameter change, so both pick up the same slope, and the flatter
#          criterion sheds it more slowly and travels further.
#
#          The criterion is a plain quadratic, not a fitted model. Nothing here
#          reads a data object or runs an optimization, so it is cheap and
#          belongs in the data tier of _run_all.R.
#
# Creator: Matthew LH. Cheng

slope_gain <- 0.3   # slope left at the old solution after a parameter changes
shed_sharp <- 2.0   # slope shed per unit of x, sharply bending criterion
shed_flat <- 0.5    # slope shed per unit of x, flat criterion

# One panel: fitted criterion, the deformed criterion, the slope left behind at
# the old solution, and the walk to the new one.
draw_panel <- function(shed, gain) {

  x <- seq(-1.6, 1.6, length.out = 800)
  h_fit <- 0.5 * shed * x^2            # criterion at the estimates
  h_new <- 0.5 * shed * x^2 + gain * x # criterion after p_j changes
  walk <- -gain / shed                 # where the new solution sits

  plot(x, h_fit, type = "l", lwd = 2.6, col = "gray20", ylim = c(-0.75, 2.2),
       xlab = expression(x == log(F)), ylab = "criterion  h", xaxt = "n")
  lines(x, h_new, lwd = 2.2, col = "gray45", lty = 2)
  abline(h = 0, col = "gray88", lwd = 0.8)

  # slope of the deformed criterion at the old solution
  x_tan <- seq(-0.55, 0.55, length.out = 2)
  lines(x_tan, gain * x_tan, lwd = 2.4, col = "steelblue4")

  points(0, 0, pch = 19, cex = 1.35, col = "gray20")
  points(walk, 0.5 * shed * walk^2 + gain * walk, pch = 19, cex = 1.35, col = "firebrick")

  arrows(0, -0.42, walk, -0.42, length = 0.07, code = 2, lwd = 2.2, col = "firebrick")
  text(walk / 2, -0.66, paste("walk", sprintf("%.2f", abs(walk))), col = "firebrick")
  axis(1, at = 0, labels = expression(x^"*"))
  mtext(bquote("slope changes" ~ .(shed) ~ "per unit of x"), side = 3, line = 0.5, cex = 0.98)

}

png("vignettes/figures/i_refpt_sensitivity.png", width = 2000, height = 900, res = 230)
par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.4, 1.0), las = 1, cex.lab = 1.05)
draw_panel(shed_sharp, slope_gain)
draw_panel(shed_flat, slope_gain)
dev.off()
