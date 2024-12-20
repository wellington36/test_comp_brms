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
