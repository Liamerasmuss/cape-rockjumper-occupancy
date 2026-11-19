# this script basically is the same as simulate_occ_two_idicators, 
# but I run the simulation 100 times to check the model actually performs well!
# corresponds to the notes written with Murray on 14 Aug (in file)


library(cmdstanr)
library(posterior)

# these change according to hypothesis #####
S <- 500
K <- 4
n_rep <- 100

psi_true <- 0.5
p_true <- 0.6
beta0_true <- -0.5
beta1_true <- 1.2
beta2_true <- 0.8

p_ind1_occ   <- 0.9
p_ind2_occ   <- 0.8
p_ind1_unocc <- 0.6
p_ind2_unocc <- 0.5
###############################################

mod <- cmdstan_model("stan/two_indicator_model.stan")

# added in medians here now to test theory about posterior median
results <- data.frame(
  rep = 1:n_rep,
  psi_mean = NA, psi_median = NA, psi_lo = NA, psi_hi = NA, psi_rhat = NA,
  p_mean = NA, p_median = NA, p_lo = NA, p_hi = NA, p_rhat = NA,
  beta0_mean = NA, beta0_median = NA, beta0_lo = NA, beta0_hi = NA, beta0_rhat = NA,
  beta1_mean = NA, beta1_median = NA, beta1_lo = NA, beta1_hi = NA, beta1_rhat = NA,
  beta2_mean = NA, beta2_median = NA, beta2_lo = NA, beta2_hi = NA, beta2_rhat = NA,
  n_divergent = NA
)

for (r in 1:n_rep) {
  
  set.seed(r)
  
  z <- rbinom(S, 1, psi_true)
  
  y1 <- matrix(NA, nrow = S, ncol = K)
  y2 <- matrix(NA, nrow = S, ncol = K)
  for (i in 1:S) {
    if (z[i] == 1) {
      p_ind1_i <- p_ind1_occ
      p_ind2_i <- p_ind2_occ
    } else {
      p_ind1_i <- p_ind1_unocc
      p_ind2_i <- p_ind2_unocc
    }
    y1[i, ] <- rbinom(K, 1, p_ind1_i)
    y2[i, ] <- rbinom(K, 1, p_ind2_i)
  }
  
  theta <- matrix(NA, nrow = S, ncol = K)
  a <- matrix(NA, nrow = S, ncol = K)
  y_crj <- matrix(NA, nrow = S, ncol = K)
  
  for (i in 1:S) {
    for (j in 1:K) {
      theta[i, j] <- plogis(beta0_true + beta1_true * y1[i, j] + beta2_true * y2[i, j])
      a[i, j] <- rbinom(1, 1, theta[i, j])
      y_crj[i, j] <- rbinom(1, 1, z[i] * a[i, j] * p_true)
    }
  }

  # for stanley
  storage.mode(y_crj) <- "integer"
  storage.mode(y1) <- "integer"
  storage.mode(y2) <- "integer"
  
  stan_data <- list(S = S, K = K, y = y_crj, y1 = y1, y2 = y2)
  
  fit <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = r,
    refresh = 0,
    show_messages = FALSE
  )
  # update: added in median here as well
  summ <- fit$summary(c("psi", "p", "beta0", "beta1", "beta2"), "mean", "median","rhat",
                      ~quantile(.x, probs = c(0.05, 0.95)))
  
  diag <- fit$diagnostic_summary()
  
  results$psi_mean[r] <- summ$mean[summ$variable == "psi"]
  results$psi_median[r] <- summ$median[summ$variable == "psi"]
  results$psi_lo[r] <- summ$`5%`[summ$variable == "psi"]
  results$psi_hi[r] <- summ$`95%`[summ$variable == "psi"]
  results$psi_rhat[r] <- summ$rhat[summ$variable == "psi"]
  
  results$p_mean[r] <- summ$mean[summ$variable == "p"]
  results$p_median[r] <- summ$median[summ$variable == "p"]
  results$p_lo[r] <- summ$`5%`[summ$variable == "p"]
  results$p_hi[r] <- summ$`95%`[summ$variable == "p"]
  results$p_rhat[r] <- summ$rhat[summ$variable == "p"]
  
  results$beta0_mean[r] <- summ$mean[summ$variable == "beta0"]
  results$beta0_median[r] <- summ$median[summ$variable == "beta0"]
  results$beta0_lo[r] <- summ$`5%`[summ$variable == "beta0"]
  results$beta0_hi[r] <- summ$`95%`[summ$variable == "beta0"]
  results$beta0_rhat[r] <- summ$rhat[summ$variable == "beta0"]
  
  results$beta1_mean[r] <- summ$mean[summ$variable == "beta1"]
  results$beta1_median[r] <- summ$median[summ$variable == "beta1"]
  results$beta1_lo[r] <- summ$`5%`[summ$variable == "beta1"]
  results$beta1_hi[r] <- summ$`95%`[summ$variable == "beta1"]
  results$beta1_rhat[r] <- summ$rhat[summ$variable == "beta1"]
  
  results$beta2_mean[r] <- summ$mean[summ$variable == "beta2"]
  results$beta2_median[r] <- summ$median[summ$variable == "beta2"]
  results$beta2_lo[r] <- summ$`5%`[summ$variable == "beta2"]
  results$beta2_hi[r] <- summ$`95%`[summ$variable == "beta2"]
  results$beta2_rhat[r] <- summ$rhat[summ$variable == "beta2"]
  
  results$n_divergent[r] <- sum(diag$num_divergent)
  
  cat("replicate", r, "of", n_rep, "done\n") # remember to always do this...super useful! (unnati, 2026)
}

results$psi_cover <- (psi_true >= results$psi_lo) & (psi_true <= results$psi_hi)
results$p_cover <- (p_true >= results$p_lo) & (p_true <= results$p_hi)
results$beta0_cover <- (beta0_true >= results$beta0_lo) & (beta0_true <= results$beta0_hi)
results$beta1_cover <- (beta1_true >= results$beta1_lo) & (beta1_true <= results$beta1_hi)
results$beta2_cover <- (beta2_true >= results$beta2_lo) & (beta2_true <= results$beta2_hi)

coverage <- c(
  psi = mean(results$psi_cover),
  p = mean(results$p_cover),
  beta0 = mean(results$beta0_cover),
  beta1 = mean(results$beta1_cover),
  beta2 = mean(results$beta2_cover)
)

print(coverage)

saveRDS(results, "output/recovery_two_indicators_S500_baseline_mark2.rds")
