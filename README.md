
bdist
=====

This is a package of distribution functions I have found useful for Bayesian analysis. Distributions include:

-   double Generalized Pareto
-   Weibull (alternative parameterization)
-   student t (alternative parameterization)
-   truncated normal
-   multivariate normal
-   inverse gamma
-   inverse root gamma

All distributions are fully unit tested.

Note: This pacakge is intended primarily for myself. Use it at your own risk!

Installation
------------

    library(devtools)
    install_github('kuperov/bdist')

Example
-------

Plot a double generalized Pareto distribution:

``` r
  library(bdist)
  gdp.density <- function(x) dgdp(x, xi = 1, alpha = 2)
  plot(gdp.density, from = -5, to = 5, col = 'blue',
       main = expression("Double Generalized Pareto distribution:" ~
                           xi == 1~alpha == 2))
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)
