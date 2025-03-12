# Log of the k-th term of the COM-Poisson normalizing constant
log_k_term <- function(log_mu, nu, k) {
  (k - 1) * log_mu - nu * lgamma(k)
}

# Bound for the remainder of the COM-Poisson normalizing constant
bound_remainder <- function(k_current_term, k_previous_term) {
  k_current_term - log(-expm1(k_current_term - k_previous_term))
}

# Log normalizing constant of the COM-Poisson distribution
log_Z_com_poisson <- function(log_mu, nu, M = 10000, leps = -52 * log(2)) {
  if (nu <= 0 || is.infinite(nu)) stop("nu must be positive and finite")
  
  log_Z_terms <- numeric(M)
  log_Z_terms[1] <- log_k_term(log_mu, nu, 1)
  log_Z_terms[2] <- log_k_term(log_mu, nu, 2)
  
  k <- 2
  while ((log_Z_terms[k] >= log_Z_terms[k - 1]) ||
         (bound_remainder(log_Z_terms[k], log_Z_terms[k - 1]) >= leps) &&
         (k < M)) {
    k <- k + 1
    log_Z_terms[k] <- log_k_term(log_mu, nu, k)
  }
  
  log_sum_exp <- function(x) max(x) + log(sum(exp(x - max(x))))
  log_Z <- log_sum_exp(log_Z_terms[1:k])
  
  # cat("k =", k, "Z =", exp(log_Z), "\n")
  return(log_Z)
}

# COM-Poisson log-likelihood
com_poisson_log_lpmf <- function(y, log_mu, nu) {
  if (nu == 1) return(dpois(y, exp(log_mu), log = TRUE))
  y * log_mu - nu * lgamma(y + 1) - log_Z_com_poisson(log_mu, nu)
}

# Log-posterior
log_posterior <- function(y, intercept, shape) {
  log_mu <- intercept
  
  # Prior: Intercept ~ Student-t(3, -2.3, 2.5) and shape ~ Gamma(0.01, 0.01)
  t_density <- dt((intercept + 2.3) / 2.5, df = 3, log = TRUE) - log(2.5)
  prior <- t_density + dgamma(shape, shape = 0.01, rate = 0.01, log = TRUE)
  
  likelihood <- sum(sapply(y, com_poisson_log_lpmf, log_mu = log_mu, nu = shape))
  
  return(prior + likelihood)
}

# Metropolis-Hastings MCMC
mcmc <- function(y, n_iter = 1000) {
  intercept <- 0
  shape <- 1
  
  chain <- matrix(NA, nrow = n_iter, ncol = 2)
  colnames(chain) <- c("intercept", "shape")
  
  for (i in 1:n_iter) {
    # Propose new values
    intercept_prop <- rnorm(1, mean = intercept, sd = 0.1)
    shape_prop <- abs(rnorm(1, mean = shape, sd = 0.1))
    
    # Calculate acceptance ratio
    log_r <- log_posterior(y, intercept_prop, shape_prop) -
      log_posterior(y, intercept, shape)
    
    if (log(runif(1)) < log_r) {
      intercept <- intercept_prop
      shape <- shape_prop
    }
    
    chain[i, ] <- c(intercept, shape)
    
    if (i %% 50 == 0) cat("Iteration", i, "\n")
  }
  
  return(chain)
}

# Load data
zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")

# Run MCMC
set.seed(666)
chain <- mcmc(zinb$count, n_iter = 1000)

# Inspect results
plot(chain[, 1], type = "l", main = "Intercept")
plot(chain[, 2], type = "l", main = "Shape")

summary(chain)
