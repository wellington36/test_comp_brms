# Test new COMPoisson of brms

Running `simple_test.R`

```r
> print(results)
        parameters COMP_updated COMP_brms COMP_true COMPoissonReg
1 mu=0.50, nu=2.00     1.266066  1.063483  1.266066      1.266066
2 mu=1.00, nu=1.50     2.430915  2.430915  2.430915      2.430915
3 mu=1.10, nu=1.40     2.781602  2.926655  2.781602      2.781601
4 mu=2.00, nu=1.30     8.193377 14.437926  8.200838      8.200835
5 mu=3.00, nu=1.20    25.058946 59.294582 25.066949     25.066940
> print(results_error)
        parameters COMP_updated COMP_poisson_brms COMPoissonReg
1 mu=0.50, nu=2.00  0.000000000         0.2025825  6.829020e-08
2 mu=1.00, nu=1.50  0.000000000         0.0000000  1.282382e-07
3 mu=1.10, nu=1.40  0.000000000         0.1450533  1.094493e-06
4 mu=2.00, nu=1.30  0.007460465         6.2370881  2.932471e-06
5 mu=3.00, nu=1.20  0.008003575        34.2276329  9.000607e-06
```

We also do a test with MCMC, first install my brms's version:

```R
devtools::install_github("wellington36/brms")
```

then MCMC can be run:

```R
> library(brms)
> zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")
>
> fit <- brm(count ~ 1,
             cores=4,
             iter=1000,
             data = zinb,
             family = "com_poisson",
             control = list(adapt_delta = 0.99))
...
> summary(fit)
 Family: com_poisson 
  Links: mu = log; shape = identity 
Formula: count ~ 1 
   Data: zinb (Number of observations: 250) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept    -0.27      0.02    -0.30    -0.24 1.00      407      319

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.00      0.00     0.00     0.00 1.01      208      373

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
Warning message:
There were 53 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
```

And with `plot(fit)`

![image](https://github.com/user-attachments/assets/677e670a-ce72-4765-8ea4-891ef857c174)
