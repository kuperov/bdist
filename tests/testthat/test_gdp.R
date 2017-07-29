context('generalized double Pareto distribution')

set.seed(1234) # can comment out for non-automated tests

params <- data.frame(
  alpha=c(1.5,2,2.5), xi=c(0.1, 1, 10))

integral.tol <- 1e-4

for (i in 1:nrow(params)) {
  alpha = params$alpha[i]
  xi = params$xi[i]

  test_that(sprintf('GDP density integrates to 1 (alpha=%f, xi=%f)', alpha, xi), {
    igral <- integrate(function(x) dgdp(x, alpha = alpha, xi = xi), lower=-Inf, upper=Inf)
    expect_lt(abs(1 - igral$value), integral.tol)
  })

  limits <- c(0.5, 1.5, 50, 100)
  lapply(limits, function(lim) {
    msg <- sprintf('GDP density integrates to cdf (alpha=%f, xi=%f, 0 to %f)', alpha, xi, lim)
    test_that(msg, {
      igral <- integrate(function(x) dgdp(x, alpha = alpha, xi = xi), lower=-Inf, upper=lim)
      p <- pgdp(lim, alpha = alpha, xi = xi)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('KS test for GDP & cdf (alpha=%f, xi=%f)', alpha, xi), {
    draws <- rgdp(1e5, alpha = alpha, xi = xi)
    kt <- ks.test(draws, 'pgdp', alpha = alpha, xi = xi)
    expect_gt(kt$p.value, 0.01)
  })
}

test_that('GDP draws the right number of variates', {
  ws <- rgdp(n = 243, alpha = 1.3, xi = 0.01)
  expect_equal(243, length(ws))
})
