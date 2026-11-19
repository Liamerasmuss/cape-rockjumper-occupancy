// Simplest single-season occupancy model (constant detection), ragged visits.
// Same structure and likelihood code as the composite-indicator model, without access.
functions {
  real lse2(real a, real b) {
    if (a == negative_infinity() && b == negative_infinity()) {
      return negative_infinity();
    }
    return log_sum_exp(a, b);
  }
}
data {
  int<lower=1> S;
  int<lower=1> Kmax;
  array[S] int<lower=1, upper=Kmax> K;
  array[S, Kmax] int<lower=0, upper=1> y;   // padded with 0 beyond K[i], never used
}
parameters {
  real<lower=0, upper=1> psi;
  real<lower=0, upper=1> p;   // detection per visit to an occupied pentad (no access term)
}
model {
  psi ~ uniform(0, 1);
  p ~ uniform(0, 1);

  for (i in 1:S) {
    real ell1 = log(psi);
    int y_sum = 0;
    for (j in 1:K[i]) {
      ell1 += bernoulli_lpmf(y[i, j] | p);
      y_sum += y[i, j];
    }
    real ell0 = (y_sum == 0) ? log1m(psi) : negative_infinity();
    target += lse2(ell0, ell1);
  }
}
