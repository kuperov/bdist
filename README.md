
acfunc
======

This is a package of distribution functions I have found useful for Bayesian analysis. It isn't really intended for public release: others are welcome to use this code, but be aware that you do so at your own risk.

Installation
------------

    library(devtools)
    install_github('kuperov/acfunc')

Example
-------

Plot a Weibull distribution

``` r
  library(acfunc)
  gdp.density <- function(x) dgdp(x, xi = 1, alpha = 2)
  plot(gdp.density, from = -5, to = 5, col = 'blue',
       main = expression("Double Generalized Pareto distribution:" ~
                           xi == 1~alpha == 2))
```

![](README_files/figure-markdown_github/unnamed-chunk-1-1.png)
