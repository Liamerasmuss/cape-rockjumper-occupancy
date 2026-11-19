functions {
  // safe pairwise log_sum_exp: if BOTH branches are -inf (impossible),
  // return -inf directly instead of letting log_sum_exp produce NaN gradients
  // Claude used for this as I kept getting an error and did not know hot to solve it
  // 
  real lse2(real a, real b) {
    if (a == negative_infinity() && b == negative_infinity()) {
      return negative_infinity();
    }
    return log_sum_exp(a, b);
  }
}

data {
  int<lower=1> S;  // number of sites
  int<lower=1> K;   // number of visits per site // vector of length s 
  array[S, K] int<lower=0, upper=1> y; // CRJ detection histories
  array[S, K] int<lower=0, upper=1> y1;   // aux species 1 detections , max no visits
  array[S, K] int<lower=0, upper=1> y2; // aux species 2 detections
}

parameters {
  real<lower=0, upper=1> psi;
  real<lower=0, upper=1> p;
  real beta0;    // intercept for theta_ij (acc probability)
  real beta1; // effect of aux species 1 on access prob
  real beta2; // effect of aux species 2 on access prob
}

// transformed parameters {
  // back-transform to the prob scale 
  // real psi = inv_logit(logit_psi);
  // real p = inv_logit(logit_p);
// }

model {
  // Murray informed these
  psi ~ uniform(0, 1);
  p ~ uniform(0, 1);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(1, 0.5); // fine
  beta2 ~ normal(1, 0.5);

  for (i in 1:S) {

    real ell1 = log(psi);
    int y_sum = 0;
    for (j in 1:K) { // for j in 1:ki
      real theta_ij = inv_logit(beta0 + beta1 * y1[i, j] + beta2 * y2[i, j]);
      real theta_p = theta_ij * p;
      ell1 += bernoulli_lpmf(y[i, j] | theta_p);
      y_sum += y[i, j];
    }
    //an unoccupied site can't produce a crj detection
    real ell0 = (y_sum == 0) ? log1m(psi) : negative_infinity();

    target += lse2(ell0, ell1);
  }
}
