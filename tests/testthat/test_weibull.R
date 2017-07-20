context('weibull distribution')

set.seed(1234) # can comment out for non-automated tests

params <- data.frame(
  alpha=c(1.5,2,2.5), gamma=c(0.01, 0.1, 1))

integral.tol <- 1e-4

for (i in 1:nrow(params)) {
  alpha = params$alpha[i]
  gamma = params$gamma[i]

  test_that(sprintf('Weibull density integrates to 1 (alpha=%f, gamma=%f)', alpha, gamma), {
    igral <- integrate(function(x) dweibull2(x, alpha = alpha, gamma = gamma), lower=0, upper=Inf)
    expect_lt(abs(1 - igral$value), integral.tol)
  })

  limits <- c(0.5, 1.5, 50, 100)
  lapply(limits, function(lim) {
    msg <- sprintf('Weibull density integrates to cdf (alpha=%f, gamma=%f, 0 to %f)', alpha, gamma, lim)
    test_that(msg, {
      igral <- integrate(function(x) dweibull2(x, alpha = alpha, gamma = gamma), lower=0, upper=lim)
      p <- pweibull2(lim, alpha = alpha, gamma = gamma)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('KS test for Weibull & cdf (alpha=%f, gamma=%f)', alpha, gamma), {
    draws <- rweibull2(1e4, alpha = alpha, gamma = gamma)
    kt <- ks.test(draws, 'pweibull2', alpha = alpha, gamma = gamma)
    expect_gt(kt$p.value, 0.01)
  })
}

test_that('Weibull draws the right number of variates', {
  ws <- rweibull2(n = 243, alpha = 1.3, gamma = 0.01)
  expect_equal(243, length(ws))
})
