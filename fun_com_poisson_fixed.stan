functions {
// log of kth term of the normalizing series of the COM Poisson distribution
real log_k_term(real log_mu, real nu, int k) {
  return (k - 1) * log_mu - nu * lgamma(k);
}

// log normalizing constant of the COM Poisson distribution
// implementation inspired by code of Ben Goodrich
// improved following suggestions of Sebastian Weber (#892)
// Args:
//   log_mu: log location parameter
//   shape: positive shape parameter
real fixed_log_Z_com_poisson(real log_mu, real nu) {
  real log_Z;
  int k = 2;
  int M = 1000000;
  vector[M] log_Z_terms;

  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  
  // first 2 terms of the series
  log_Z_terms[1] = log_k_term(log_mu, nu, 1);
  log_Z_terms[2] = log_k_term(log_mu, nu, 2);

  while (k < M) {
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
real fixed_com_poisson_log_lpmf(int y, real log_mu, real nu) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  return y * log_mu - nu*lgamma(y + 1) - fixed_log_Z_com_poisson(log_mu, nu);
}
// COM Poisson log-PMF for a single response
real fixed_com_poisson_lpmf(int y, real mu, real nu) {
  if (nu == 1) return poisson_lpmf(y | mu);
  return fixed_com_poisson_log_lpmf(y | log(mu), nu);
}
// COM Poisson log-CDF for a single response
real fixed_com_poisson_lcdf(int y, real mu, real nu) {
  real log_mu;
  real log_Z;  // log denominator
  vector[y] log_num_terms;  // terms of the log numerator
  if (nu == 1) {
    return poisson_lcdf(y | mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  if (y > 10000) {
    reject("cannot handle y > 10000");
  }
  log_mu = log(mu);
  if (y * log_mu - nu*lgamma(y + 1) <= -36.0) {
    // y is large enough for the CDF to be very close to 1;
    return 0;
  }
  log_Z = fixed_log_Z_com_poisson(log_mu, nu);
  if (y == 0) {
    return -log_Z; 
  }
  // first 2 terms of the series
  log_num_terms[1] = log1p_exp(nu * log_mu);
  // remaining terms of the series until y
  for (k in 2:y) {
    log_num_terms[k] = k * log_mu - nu*lgamma(k + 1);
  }
  return log_sum_exp(log_num_terms) - log_Z;
}
// COM Poisson log-CCDF for a single response
real fixed_com_poisson_lccdf(int y, real mu, real nu) {
  return log1m_exp(fixed_com_poisson_lcdf(y | mu, nu));   
}
}