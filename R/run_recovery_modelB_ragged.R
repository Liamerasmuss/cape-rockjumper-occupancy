# Recovery study for Model B: composite indicator + ragged visits.
# Data are simulated from Model B's generative process, using the REAL number of
# visits per pentad (so the simulated data is ragged in the same way), then Model B is
# fitted to each simulated dataset and the truth is compared with the estimates.
#
# indicator_ratio controls how the indicator relates to true occupancy:
#   1 -> indicator independent of occupancy (matches the model exactly: tests the CODE)
#   3 -> indicator 3x more likely at occupied sites (co-occurrence: tests ROBUSTNESS)
#
# Judge the study mainly on psi, q0 and q1, where q_s = plogis(beta0 + beta1 * s) * p.
# p, beta0 and beta1 are not separately identified, so do not expect clean recovery.
#
# Run from the project root (the folder containing data/, stan/ and output/).

library(cmdstanr)

# ---- settings ----

n_rep <- 100              # try 2 first to check the script runs, then raise
indicator_ratio <- 1     # 1 = independent of occupancy, 3 = co-occurrence

psi_true <- 0.25
p_true <- 0.92
beta0_true <- -2.5
beta1_true <- 3.0

s_rate <- 0.0242         # overall rate of the composite indicator in the real data

# implied indicator probability at unoccupied / occupied sites, keeping the overall
# rate at s_rate:  s_rate = psi * p_occ + (1 - psi) * p_unocc,  p_occ = ratio * p_unocc
p_s_unocc <- s_rate / (psi_true * indicator_ratio + (1 - psi_true))
p_s_occ <- indicator_ratio * p_s_unocc

q0_true <- plogis(beta0_true) * p_true
q1_true <- plogis(beta0_true + beta1_true) * p_true

pars <- c("psi", "p", "beta0", "beta1", "q0", "q1")
truth <- c(psi = psi_true, p = p_true, beta0 = beta0_true, beta1 = beta1_true,
           q0 = q0_true, q1 = q1_true)

# ---- real visit counts per pentad ----

raw <- readRDS("data/clean/cleaned_multispp_occ_data.RDS")  # adjust filename if different
K <- as.integer(table(raw$Pentad))
S <- length(K)
Kmax <- max(K)
cat("S =", S, " Kmax =", Kmax, " median K =", median(K), "\n")

mod <- cmdstan_model("stan/composite_indicator_model_ragged.stan")

# ---- storage ----

est_mean <- matrix(NA, nrow = n_rep, ncol = length(pars), dimnames = list(NULL, pars))
est_lo <- est_mean
est_hi <- est_mean
max_rhat <- rep(NA, n_rep)
n_divergent <- rep(NA, n_rep)

# ---- replicates ----

for (r in 1:n_rep) {

  set.seed(r)

  y <- matrix(0L, nrow = S, ncol = Kmax)
  synth <- matrix(0L, nrow = S, ncol = Kmax)

  for (i in 1:S) {
    k_i <- K[i]
    z_i <- rbinom(1, 1, psi_true)

    p_s <- if (z_i == 1) p_s_occ else p_s_unocc
    s_i <- rbinom(k_i, 1, p_s)

    theta_i <- plogis(beta0_true + beta1_true * s_i)
    a_i <- rbinom(k_i, 1, theta_i)
    y_i <- rbinom(k_i, 1, z_i * a_i * p_true)

    synth[i, 1:k_i] <- s_i
    y[i, 1:k_i] <- y_i
  }

  if (r == 1) {
    cat("Simulated indicator rate:", round(sum(synth) / sum(K), 4), "\n")
    cat("Simulated detection rate:", round(sum(y) / sum(K), 4), "\n")
  }

  stan_data <- list(S = S, Kmax = Kmax, K = K, y = y, synth = synth)

  fit <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    seed = r,
    refresh = 0,
    show_messages = FALSE
  )

  draws <- as.data.frame(fit$draws(c("psi", "p", "beta0", "beta1"), format = "df"))
  draws$q0 <- plogis(draws$beta0) * draws$p
  draws$q1 <- plogis(draws$beta0 + draws$beta1) * draws$p

  for (nm in pars) {
    est_mean[r, nm] <- mean(draws[[nm]])
    qs <- quantile(draws[[nm]], probs = c(0.05, 0.95))
    est_lo[r, nm] <- qs[1]
    est_hi[r, nm] <- qs[2]
  }

  summ <- fit$summary(c("psi", "p", "beta0", "beta1"), "rhat")
  max_rhat[r] <- max(summ$rhat)

  diag <- fit$diagnostic_summary()
  n_divergent[r] <- sum(diag$num_divergent)

  cat("replicate", r, "of", n_rep, "done\n")
}

# ---- summaries ----

coverage <- sapply(pars, function(nm) mean(truth[nm] >= est_lo[, nm] & truth[nm] <= est_hi[, nm]))
bias <- sapply(pars, function(nm) mean(est_mean[, nm] - truth[nm]))
width <- sapply(pars, function(nm) mean(est_hi[, nm] - est_lo[, nm]))

print(round(data.frame(truth = truth, coverage = coverage, bias = bias, width = width), 3))
cat("\nReplicates with divergences:", sum(n_divergent > 0), "\n")
cat("Replicates with max R-hat > 1.01:", sum(max_rhat > 1.01), "\n")

saveRDS(
  list(truth = truth, est_mean = est_mean, est_lo = est_lo, est_hi = est_hi,
       max_rhat = max_rhat, n_divergent = n_divergent, indicator_ratio = indicator_ratio),
  file = sprintf("output/recovery_100_modelB_ragged_ratio%s.rds", indicator_ratio)
)
max(max_rhat)
