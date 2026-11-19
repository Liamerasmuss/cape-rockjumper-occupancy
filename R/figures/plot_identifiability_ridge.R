# Identifiability ridge figure for Model B (original and wide prior side by side).
# Joint posterior of p and exp(beta0), with the curve of constant q0 = theta0 * p
# drawn over it. The data pin down q0, so the posterior lies along this curve.
# Uses the saved draws, no refitting. Run from the project root.

orig <- readRDS("output/realfit_modelB_normal_draws.rds")
wide <- readRDS("output/realfit_modelB_wide_draws.rds")

# q0 = plogis(beta0) * p = p * x / (1 + x) with x = exp(beta0), so p = q0 * (1 + x) / x
ridge <- function(q0, x) q0 * (1 + x) / x                      # may exceed 1
ridge_clip <- function(q0, x) pmin(q0 * (1 + x) / x, 1)       # for the shaded band

panel <- function(d, title, xmin, xmax) {
  x <- exp(d$beta0)
  q0_lo <- quantile(d$q0, 0.05)
  q0_hi <- quantile(d$q0, 0.95)
  q0_md <- median(d$q0)
  r <- cor(d$p, x)
  
  xs <- seq(xmin, xmax, length.out = 400)
  plot(NA, xlim = c(xmin, xmax), ylim = c(0.4, 1.02), xlab = expression(e^beta[0]),
       ylab = expression(p), bty = "l", main = title, font.main = 1)
  
  # band of q0 values with 90% posterior mass
  polygon(c(xs, rev(xs)), c(ridge_clip(q0_lo, xs), rev(ridge_clip(q0_hi, xs))),
          col = adjustcolor("firebrick", 0.18), border = NA)
  
  # posterior draws (thinned so the PDF stays small)
  idx <- seq(1, nrow(d), length.out = min(nrow(d), 3000))
  points(x[idx], d$p[idx], pch = 16, cex = 0.45, col = adjustcolor("grey20", 0.25))
  
  r_md <- ridge(q0_md, xs); r_md[r_md > 1] <- NA
  lines(xs, r_md, col = "firebrick", lwd = 2)
  abline(h = 1, lty = 3, col = "grey50")
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c("Posterior draws",
                    expression(q[0] == theta[0] * p ~ "(posterior median)"),
                    expression("90% interval for" ~ q[0])),
         pch = c(16, NA, 15), lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
         col = c(adjustcolor("grey20", 0.6), "firebrick", adjustcolor("firebrick", 0.3)),
         pt.cex = c(0.8, NA, 1.4))
  mtext(sprintf("cor = %.2f", r), side = 3, line = -1.3, adj = 0.97, cex = 0.85)
}

# shared x-range, chosen to cover both posteriors (check with range(exp(orig$beta0)), range(exp(wide$beta0)))
xmin <- 0.05
xmax <- 0.17
dir.create("figs", showWarnings = FALSE)
png("figs/identifiability-ridge.png", width = 2200, height = 1000, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), mgp = c(2.4, 0.7, 0), las = 1)
panel(orig, expression("Original prior," ~ beta[1] ~ "~ N(1, " * 0.5^2 * ")"), xmin, xmax)
panel(wide, expression("Wide prior," ~ beta[1] ~ "~ N(0, " * 2^2 * ")"), xmin, xmax)
dev.off()