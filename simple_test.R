library(rstan)
library(COMPoissonReg)

Rcpp::cppFunction('double logDiffExp(double x, double y)
{return x > y ? 
Rf_logspace_sub(x, y) :
Rf_logspace_sub(y, x);}')

expose_stan_functions("fun_com_poisson_updated.stan")
expose_stan_functions("fun_com_poisson_brms.stan")
expose_stan_functions("fun_com_poisson_fixed.stan")

# considering the parametrization of https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution
y <- 0:10
mu <- 1.1
nu <- 1.1
lambda <- mu^nu


updated <- numeric(length(y))
brms <- numeric(length(y))
fixed <- numeric(length(y))

log_Z_updated <- log_Z_com_poisson(log(lambda), nu)
log_Z_brms <- brms_log_Z_com_poisson(log(lambda), nu)
log_Z_fixed <- fixed_log_Z_com_poisson(log(lambda), nu)

for (i in seq_along(y)) {
  updated[i] <- exp(com_poisson_lpmf(y[i], lambda, nu))
}

for (i in seq_along(y)) {
  brms[i] <- exp(brms_com_poisson_lpmf(y[i], lambda, nu))
}

for (i in seq_along(y)) {
  fixed[i] <- exp(fixed_com_poisson_lpmf(y[i], lambda, nu))
}

results <- data.frame(
  y = y,
  com_poisson_updated = updated,
  com_poisson_brms = brms,
  com_poisson_fixed = fixed,
  COMPoissonReg = dcmp(y, mu, nu)
)

print(results)
print(exp(logDiffExp(log_Z_updated, log_Z_fixed)))
print(exp(logDiffExp(log_Z_brms, log_Z_fixed)))