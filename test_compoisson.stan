functions {
real log_Z_com_poisson(real log_mu, real nu) {
  
}

real com_poisson_log_lpmf(int y, real log_mu, real nu) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  return y * log_mu - nu*lgamma(y + 1) - log_Z_com_poisson(log_mu, nu);
}

real com_poisson_lpmf(int y, real mu, real nu) {
  if (nu == 1) return poisson_lpmf(y | mu);
  return com_poisson_log_lpmf(y | log(mu), nu);
}
}


data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, -2.3, 2.5);
  lprior += gamma_lpdf(shape | 0.01, 0.01);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
      target += com_poisson_log_lpmf(Y[n] | mu[n], shape);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}