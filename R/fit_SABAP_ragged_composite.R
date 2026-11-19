
# Script for using the cleaned SABAP data and a composite indicator species 
# using ragged arrays
# then it fits composite_indicator_model_ragged.stan.

# synth_ij = 1 if (Sentinel Rock-thrush OR Ground Woodpecker seen)
#              AND (Orange-breasted Sunbird OR Cape Siskin seen)

library(cmdstanr)
library(posterior)

# ---- load cleaned data ----

raw <- readRDS("data/clean/cleaned_multispp_occ_data.RDS") 

# ---- build the synthetic / composite indicator ----

rocky <- (raw$Spp.sert != "-") | (raw$Spp.grwo != "-")
fynbos <- (raw$Spp.obsu != "-") | (raw$Spp.casi != "-")
raw$synth <- as.integer(rocky & fynbos)

# ---- build per site visit lists ----

pentads <- unique(raw$Pentad)
S <- length(pentads)

site_data <- vector("list", S)
K <- integer(S)

for (i in 1:S) {
  site_rows <- raw[raw$Pentad == pentads[i], ]
  site_rows <- site_rows[order(site_rows$StartDate), ]
  site_data[[i]] <- site_rows
  K[i] <- nrow(site_rows)
  cat("site", i, "of", S, "-- K =", K[i], "\n")
}

Kmax <- max(K)
cat("\nS =", S, " Kmax =", Kmax, "\n")
cat("Median K:", median(K), " Mean K:", round(mean(K), 1), "\n")

# ---- build padded matrices ----

y_crj <- matrix(0L, nrow = S, ncol = Kmax)
synth <- matrix(0L, nrow = S, ncol = Kmax)

for (i in 1:S) {
  rows <- site_data[[i]]
  k_i <- K[i]
  y_crj[i, 1:k_i] <- as.integer(rows$Spp.caro != "-")
  synth[i, 1:k_i] <- rows$synth
}

storage.mode(y_crj) <- "integer"
storage.mode(synth) <- "integer"

# ---- fit ----

stan_data <- list(S = S, Kmax = Kmax, K = K, y = y_crj, synth = synth)

mod <- cmdstan_model("stan/composite_indicator_model_ragged.stan")

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99, # made this to 0.99 to fix 11 divergent runs
  seed = 1
)
# murray says that if adapt_delta solved the issue, then the issue was not that big to begin with 

fit$diagnostic_summary()

summ <- fit$summary(c("psi", "p", "beta0", "beta1"))
print(summ)




draws <- as_draws_df(fit$draws(c("p", "beta0")))

odds_beta0 <- exp(draws$beta0)
p_draws <- draws$p

cor(p_draws, odds_beta0)

plot(odds_beta0, p_draws,
     xlab = "e^beta0  (baseline access odds)",
     ylab = "p  (detection probability)",
     main = "Joint posterior: p vs e^beta0",
     pch = 16, col = rgb(0, 0, 0, 0.15))

# overlay a 2D density contour on top of the scatter, to see the shape more clearly
dens <- MASS::kde2d(odds_beta0, p_draws, n = 100)
contour(dens, add = TRUE, col = "red", lwd = 1.5)
# murray says that this does not look tooo bad. YAYYYYY


#########################################################################################################
# check on thingie 
# Prior versus posterior densities for Model B (original priors).
# Assumes `fit` is the real-data fit of synthetic_indicator_model_ragged.stan
# (fitted with adapt_delta = 0.99). Run from the project root.
# For the wide-prior fit, use that fit object, change the beta priors below
# to dnorm(x, 0, 2), and change the output file name.

library(posterior)

draws <- as_draws_df(fit$draws(c("psi", "p", "beta0", "beta1")))

plot_prior_post <- function(post, prior_density, xlim, main) {
  x <- seq(xlim[1], xlim[2], length.out = 500)
  d <- density(post)
  ymax <- max(d$y, prior_density(x))
  plot(d, xlim = xlim, ylim = c(0, ymax), main = main, xlab = "", lwd = 2)
  lines(x, prior_density(x), col = "red", lwd = 2, lty = 2)
}

#png("figs/prior-posterior-wide-priors.png", width = 1600, height = 1200, res = 200)

par(mfrow = c(2, 2))

plot_prior_post(draws$psi, function(x) dunif(x, 0, 1), c(0, 1), "psi")
plot_prior_post(draws$p, function(x) dunif(x, 0, 1), c(0, 1), "p")
plot_prior_post(draws$beta0, function(x) dnorm(x, 0, 1), c(-4, 4), "beta0")
plot_prior_post(draws$beta1, function(x) dnorm(x, 1, 0.5), c(-1, 5), "beta1")
legend("topleft", legend = c("posterior", "prior"),
       col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")

dev.off()

######################################################################################################

draws <- as.data.frame(as_draws_df(fit$draws(c("psi", "p", "beta0", "beta1"))))
draws$q0 <- plogis(draws$beta0) * draws$p
draws$q1 <- plogis(draws$beta0 + draws$beta1) * draws$p
draws$ratio <- draws$q1 / draws$q0
draws$theta0 <- plogis(draws$beta0)
draws$theta1 <- plogis(draws$beta0 + draws$beta1)

tab <- t(sapply(draws[, c("psi", "p", "beta0", "beta1", "q0", "q1", "ratio", "theta0", "theta1")],
                function(x) c(mean = mean(x), median = median(x),
                              lo = unname(quantile(x, 0.05)),
                              hi = unname(quantile(x, 0.95)))))
print(round(tab, 3))

saveRDS(draws, "output/realfit_modelB_normal_draws.rds")   # 
