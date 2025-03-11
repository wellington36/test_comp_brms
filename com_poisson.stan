// generated with brms 2.22.8
functions {
// log approximate normalizing constant of the COM poisson distribuion
// based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
// Args: see log_Z_com_poisson()
real log_Z_com_poisson_approx(real log_mu, real nu) {
  real nu2 = nu^2;
  real log_common = log(nu) + log_mu/nu;
  array[4] real resids;
  real ans;
  real lcte = (nu * exp(log_mu/nu)) -
    ( (nu-1)/(2*nu)* log_mu + (nu-1)/2*log(2*pi()) + 0.5 *log(nu));
  real c_1 = (nu2-1)/24;
  real c_2 = (nu2-1)/1152*(nu2 + 23);
  real c_3 = (nu2-1)/414720* (5*square(nu2) - 298*nu2 + 11237);
  resids[1] = 1;
  resids[2] = c_1 * exp(-1 * log_common);
  resids[3] = c_2 * exp(-2 * log_common);
  resids[4] = c_3 * exp(-3 * log_common);
  ans = lcte + log(sum(resids));
  return ans;
}

// log of kth term of the normalizing series of the COM Poisson distribution
// Args:
//   log_mu: log location parameter
//   shape: positive shape parameter
//   k: k-th term
real log_k_term(real log_mu, real nu, int k) {
  return (k - 1) * log_mu - nu * lgamma(k);
}

// bound for the remainder of the normalizing series of the COM Poisson
// distribution given the last two terms in log-scale
// Args:
//   k_current_term: the log of a_k term
//   k_previous_term: the log of a_(k-1) term
real bound_remainder(real k_current_term, real k_previous_term) {
  return k_current_term - log(- expm1(k_current_term - k_previous_term));
}

// log normalizing constant of the COM Poisson distribution
// implementation inspired by code of Ben Goodrich
// improved following suggestions of Sebastian Weber (#892)
// Args:
//   log_mu: log location parameter
//   shape: positive shape parameter
real log_Z_com_poisson(real log_mu, real nu) {
  real log_Z;
  int k = 2;
  int M = 10000;
  real leps = -8 * log2();
  vector[M] log_Z_terms;

  if (nu == 1) {
    return exp(log_mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  //if (log_mu * nu >= log(1.5) && log_mu >= log(1.5)) {
  //  return log_Z_com_poisson_approx(log_mu, nu);
  //}
  // direct computation of the truncated series
  // check if the Mth term of the series pass in the stopping criteria
  if (bound_remainder(log_k_term(log_mu, nu, M),
                      log_k_term(log_mu, nu, M-1)) >= leps) {
    reject("nu is too close to zero.");
  }

  // first 2 terms of the series
  log_Z_terms[1] = log_k_term(log_mu, nu, 1);
  log_Z_terms[2] = log_k_term(log_mu, nu, 2);

  while ((log_Z_terms[k] >= log_Z_terms[k-1]) ||
    (bound_remainder(log_Z_terms[k], log_Z_terms[k-1]) >= leps) &&
    k < M) {
    k += 1;
    log_Z_terms[k] = log_k_term(log_mu, nu, k);
  }
  log_Z = log_sum_exp(log_Z_terms[1:k]);

  return log_Z;
}
// COM Poisson log-PMF for a single response (log parameterization)
// Args:
//   y: the response value
//   log_mu: log location parameter
//   shape: positive shape parameter
real com_poisson_log_lpmf(int y, real log_mu, real nu) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  print("com_poisson_log_lpmf = ", y * log_mu - nu*lgamma(y + 1) - log_Z_com_poisson(log_mu, nu));
  return y * log_mu - nu*lgamma(y + 1) - log_Z_com_poisson(log_mu, nu);
}
// COM Poisson log-PMF for a single response
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