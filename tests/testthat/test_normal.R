# test normal family of distributions

library(cubature)

set.seed(1234) # can comment out for non-automated tests

integral.tol <- 1e-4

# ---- 2D multivariate normal ----

test_that('2D multivariate normal integrates to unity', {
  # we'll use 50 as 'infinity'
  mu <- c(1,2)
  Sigma <- matrix(c(3,0.5,0.5,6), 2)
  f <- function(x)
    matrix(mvdnorm(t(x), mu=mu, Sigma=Sigma), ncol=1)
  volume <- pcubature(f, lowerLimit=c(-50, -50), upperLimit=c(50, 50), vectorInterface = TRUE)
  expect_true(abs(volume$integral - 1) < integral.tol)
})

# ---- truncated normal ----

params <- data.frame(
  mu=c(1.5,2,2.5),
  sigma=c(0.1, 1, 10),
  lower=c(-Inf, -5, 0),
  upper=c(Inf, Inf, 5))

for (i in 1:nrow(params)) {
  mu = params$mu[i]
  sigma = params$sigma[i]
  lower = params$lower[i]
  upper = params$upper[i]

  test_that(sprintf('rnorm_trunc density integrates to 1 (mu=%.1f, sigma=%.1f, lower=%.1f, upper=%.1f)',
                    mu, sigma, lower, upper), {
    dens <- function(x) dnorm_trunc(x, mean = mu, sd = sigma)
    igral <- integrate(dens, lower=-Inf, upper=Inf)
    expect_lt(abs(1 - igral$value), integral.tol)
  })

  limits <- c(0.5, 1.5, 5, 10)
  lapply(limits, function(lim) {
    msg <- sprintf('dnorm_trunc density integrates to cdf (mu=%.1f, sigma=%.2f, lower=%.1f, upper=%.1f, 0 to %.1f)',
                    mu, sigma, lower, upper, lim)
    test_that(msg, {
      f <- function(x) dnorm_trunc(x, mean = mu, sd = sigma, lower = lower,
                                    upper = upper)
      igral <- integrate(f, lower = -lim, upper = lim)
      p <- pnorm_trunc(lim, mean = mu, sd = sigma, lower = lower, upper = upper) -
        pnorm_trunc(-lim, mean = mu, sd = sigma, lower = lower, upper = upper)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('dnorm_trunc density and log density agree (mu=%.1f, sigma=%.1f, lower=%.1f, upper=%.1f)',
            mu, sigma, lower, upper), {
      # space out along support, cut off at 1e4 for infinite cases
      xs <- seq(max(-1e4, lower), min(1e4, upper), length.out = 10)
      dens <- dnorm_trunc(xs, mean = mu, sd = sigma, log = FALSE)
      log.dens <- dnorm_trunc(xs, mean = mu, sd = sigma, log = TRUE)
      expect_equal(dens, exp(log.dens))
    })

  test_that(sprintf('KS test for rnorm_trunc & cdf (mu=%.1f, sigma=%.1f, lower=%.1f, upper=%.1f)',
                    mu, sigma, lower, upper), {
    draws <- rnorm_trunc(1e4, mean = mu, sd = sigma, lower = lower, upper = upper)
    kt <- ks.test(draws, 'pnorm_trunc', mean = mu, sd = sigma, lower = lower, upper = upper)
    expect_gt(kt$p.value, 0.01)
  })
}

test_that('rnorm_trunc draws the right number of variates', {
  ws <- rnorm_trunc(n = 243, mean = 1.3, sd = 0.01, lower = -3, upper = 6)
  expect_equal(243, length(ws))
})
