# Fit the simplest single-season occupancy model (constant p, no access) to the
# SABAP2 data, using the same pentads, visits and sampler settings as Model B.
# Murray's suggestion !

library(cmdstanr)
library(posterior)

raw <- readRDS("data/clean/cleaned_multispp_occ_data.RDS")

# ---- build the ragged detection matrix ----

by_pentad <- split(raw, raw$Pentad)
S <- length(by_pentad)
K <- as.integer(sapply(by_pentad, nrow))
Kmax <- max(K)

y <- matrix(0L, nrow = S, ncol = Kmax)
for (i in 1:S) {
  y[i, 1:K[i]] <- as.integer(by_pentad[[i]]$Spp.caro != "-")
}

cat("S =", S, " total visits =", sum(K), " Kmax =", Kmax, "\n")
cat("Raw reporting rate:", round(sum(y) / sum(K), 4), "\n")
cat("Pentads with at least one record:", sum(rowSums(y) > 0),
    " naive occupancy:", round(mean(rowSums(y) > 0), 3), "\n")

# ---- fit ----
mod <- cmdstan_model("stan/standard_occupancy_ragged.stan")

fit <- mod$sample(
  data = list(S = S, Kmax = Kmax, K = K, y = y),
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000,
  seed = 1, refresh = 0
)

fit$diagnostic_summary()
print(fit$summary(c("psi", "p"), "mean", ~quantile(.x, c(0.05, 0.95)), "rhat", "ess_bulk"))

saveRDS(as_draws_df(fit$draws(c("psi", "p"))), "output/realfit_standard_draws.rds")
