functions {
  // leaving this here because why not. Do we need it? Not sure
  real lse2(real a, real b) {
    if (a == negative_infinity() && b == negative_infinity()) {
      return negative_infinity();
    }
    return log_sum_exp(a, b);
  }
}
data {
  int<lower=1> S; // number of sites
  int<lower=1> Kmax;  // longest visit count across all sites
  array[S] int<lower=1, upper=Kmax> K;     // number of actual visits at each site
  array[S, Kmax] int<lower=0, upper=1> y; // CRJ detections
  array[S, Kmax] int<lower=0, upper=1> synth; // synthetic / composite indicator
}
parameters {
  real<lower=0, upper=1> psi;
  real<lower=0, upper=1> p;
  real beta0;
  real beta1; // effect of composite indicator on access
}
model { // staying the same 
  psi ~ uniform(0, 1);
  p ~ uniform(0, 1);
  beta0 ~ normal(0, 1);  
  beta1 ~ normal(1, 0.5);  // wide = 0, 2

  for (i in 1:S) {
    real ell1 = log(psi);
    int y_sum = 0;
    for (j in 1:K[i]) { // this is the new part 
      real theta_ij = inv_logit(beta0 + beta1 * synth[i, j]);
      real theta_p = theta_ij * p;
      ell1 += bernoulli_lpmf(y[i, j] | theta_p);
      y_sum += y[i, j];
    }
    real ell0 = (y_sum == 0) ? log1m(psi) : negative_infinity();
    target += lse2(ell0, ell1);
  }
}
