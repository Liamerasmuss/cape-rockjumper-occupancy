# script to check the parameter recovery of the 100 simulations of the two aux species model
# plots the 90% CI of the posterior mean for psi, p, beta0, beta1, and beta2
# with the true value shown as a red line
# coverage = does the red line cross teh grey?
# bias = doe steh red line cross near the middle?
# precision = how long are the grey lines?

results <- readRDS("output/recovery_two_indicators_S500.rds") # rds saved form run_recovery_two_indiactors_sim100

# vary depending on the chosen hypothesis 
psi_true <- 0.5
p_true <- 0.6
beta0_true <- -0.5
beta1_true <- 1.2
beta2_true <- 0.8

bias <- c(
  psi = mean(results$psi_mean - psi_true),
  p = mean(results$p_mean - p_true),
  beta0 = mean(results$beta0_mean - beta0_true),
  beta1 = mean(results$beta1_mean - beta1_true),
  beta2 = mean(results$beta2_mean - beta2_true)
)
# median added in to check one hypothesis
bias_median <- c(
  psi = mean(results$psi_median - psi_true),
  p = mean(results$p_median - p_true),
  beta0 = mean(results$beta0_median - beta0_true),
  beta1 = mean(results$beta1_median - beta1_true),
  beta2 = mean(results$beta2_median - beta2_true)
)

width <- c(
  psi = mean(results$psi_hi - results$psi_lo),
  p = mean(results$p_hi - results$p_lo),
  beta0 = mean(results$beta0_hi - results$beta0_lo),
  beta1 = mean(results$beta1_hi - results$beta1_lo),
  beta2 = mean(results$beta2_hi - results$beta2_lo)
)

n_bad_rhat <- sum(
  results$psi_rhat > 1.01 | results$p_rhat > 1.01 | results$beta0_rhat > 1.01 |
    results$beta1_rhat > 1.01 | results$beta2_rhat > 1.01
)

n_divergent_reps <- sum(results$n_divergent > 0)

print(bias)
print(bias_median)
print(width)
print(n_bad_rhat)
print(n_divergent_reps)

png("figs/recovery_two_indicators_S500_baselline_mark2.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

# should possibly add in titles and stuff at a later stage
plot_recovery <- function(mean_col, lo_col, hi_col, true_val, name) {
  ord <- order(results[[mean_col]])
  plot(1:nrow(results), results[[mean_col]][ord],
       ylim = range(c(results[[lo_col]], results[[hi_col]])),
       pch = 16, cex = 0.6, xlab = "replicate (sorted)", ylab = name,
       main = name)
  arrows(1:nrow(results), results[[lo_col]][ord],
         1:nrow(results), results[[hi_col]][ord],
         length = 0, col = "grey")
  abline(h = true_val, col = "red", lwd = 2)
}

plot_recovery("psi_mean", "psi_lo", "psi_hi", psi_true, "psi")
plot_recovery("p_mean", "p_lo", "p_hi", p_true, "p")
plot_recovery("beta0_mean", "beta0_lo", "beta0_hi", beta0_true, "beta0")
plot_recovery("beta1_mean", "beta1_lo", "beta1_hi", beta1_true, "beta1")
plot_recovery("beta2_mean", "beta2_lo", "beta2_hi", beta2_true, "beta2")

dev.off()
