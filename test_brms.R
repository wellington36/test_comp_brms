################# MAIN #####################

# library(brms)
# 
# zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")
# 
# fit <- brm(count ~ 1,
#            cores=4,
#            iter=1000,
#            data = zinb,
#            family = "poisson")
# 
# summary(fit)


################# TEST #####################

library(rstan)
library(brms)

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")

stan_data <- list(
  N = nrow(zinb),
  Y = zinb$count,
  prior_only = 0
)

fit <- stan(
  file = "com_poisson.stan",
  data = stan_data,
  iter = 200,
  cores = 4
)

print(fit)
plot(fit)

################# TEST W/POISSON ###########

# library(rstan)
# library(brms)
# 
# zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")
# 
# stan_data <- list(
#   N = nrow(zinb),
#   Y = zinb$count,
#   prior_only = 0
# )
# 
# fit <- stan(
#   file = "poisson_version_of_test_brms.stan",
#   data = stan_data,
#   iter = 1000,
#   cores = 4
# )
# 
# print(fit)
# plot(fit)