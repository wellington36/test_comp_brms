library(rstan)
library(COMPoissonReg)

Rcpp::cppFunction('double logDiffExp(double x, double y)
{return x > y ? 
Rf_logspace_sub(x, y) :
Rf_logspace_sub(y, x);}'
)

expose_stan_functions("fun_com_poisson_updated.stan")
expose_stan_functions("fun_com_poisson_brms.stan")
expose_stan_functions("fun_com_poisson_fixed.stan")

# considering the parametrization of https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution
#mu <- c(0.1, 0.2, 0.4, 1, 2.5, 5, 10)
#nu <- c(10, 5, 2.5, 1, 0.4, 0.2, 0.1)
mu <- c(0.5, 1, 1.1, 2, 3)
nu <- c(2, 1.5, 1.4, 1.3, 1.2)
lambda <- mu^nu


updated <- numeric(length(mu))
brms <- numeric(length(mu))
fixed <- numeric(length(mu))
compoissonreg_lib <- numeric(length(mu))

updated_error <- numeric(length(mu))
brms_error <- numeric(length(mu))
compoissonreg_lib_error <- numeric(length(mu))


for (i in seq_along(mu)) {
  fixed[i] <- exp(fixed_log_Z_com_poisson(log(lambda[i]), nu[i]))
}

for (i in seq_along(mu)) {
  updated[i] <- exp(log_Z_com_poisson(log(lambda[i]), nu[i]))
  updated_error[i] <- abs(updated[i] - fixed[i])
}

for (i in seq_along(mu)) {
  brms[i] <- exp(brms_log_Z_com_poisson(log(lambda[i]), nu[i]))
  brms_error[i] <- abs(brms[i] - fixed[i])
}

for (i in seq_along(mu)) {
  compoissonreg_lib[i] <- 1/dcmp(0, lambda[i], nu[i])
  compoissonreg_lib_error[i] <- abs(compoissonreg_lib[i] - fixed[i])
}

parameters = numeric(length(mu))
for (i in seq_along(mu)) {
  parameters[i] <- sprintf("mu=%.2f, nu=%.2f", mu[i], nu[i])
}

results <- data.frame(
  parameters = parameters,
  COMP_updated = updated,
  COMP_brms = brms,
  COMP_true = fixed,
  COMPoissonReg = compoissonreg_lib
)


results_error <- data.frame(
  parameters = parameters,
  COMP_updated = updated_error,
  COMP_poisson_brms = brms_error,
  COMPoissonReg = compoissonreg_lib_error
)

print(results)
print(results_error)