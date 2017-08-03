context('student t distributions')

set.seed(1234) # can comment out for non-automated tests

params <- data.frame(
  mu=c(1,-6,2,-3,32,-2,40), sig.sq=c(2,3,6,1,2,3,4), df=c(2,7,4,2,1,25,35))

integral.tol <- 1e-4

for (i in 1:nrow(params)) {
  mu = params$mu[i]
  df = params$df[i]
  sig.sq = params$sig.sq[i]

  test_that(sprintf('Student t density integrates to 1 (nu=%f, mu=%f, sigma^2=%f)', df, mu, sig.sq), {
    igral <- integrate(function(x) dst(x, df = df, location = mu, scale = sig.sq), lower=-Inf, upper=Inf)
    expect_equal(1, igral$value)
  })

  test_that(sprintf('Student t density and log density agree (nu=%f, mu=%f, sigma^2=%f)', df, mu, sig.sq), {
    xs <- seq(-100*sig.sq+mu, 100*sig.sq+mu, length.out = 10)
    dens <- dst(xs, df = df, location = mu, scale = sig.sq)
    log.dens <- dst(xs, df = df, location = mu, scale = sig.sq, log = TRUE)
    expect_equal(dens, exp(log.dens))
  })

  test_that(sprintf('Student t density integrates to cdf (nu=%f, mu=%f, sigma^2=%f)', df, mu, sig.sq), {
    limits <- seq(mu-2*sig.sq, mu+2*sig.sq, length.out = 5)
    lapply(limits, function(lim) {
      igral <- integrate(function(x) dst(x, df = df, location = mu, scale = sig.sq), lower=-Inf, upper=lim)
      p <- pst(lim, df = df, location = mu, scale = sig.sq)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('KS test for student t & cdf (nu=%f, mu=%f, sigma^2=%f)', df, mu, sig.sq), {
    draws <- rst(1e4, df = df, location = mu, scale = sig.sq)
    kt <- ks.test(draws, 'pst', df = df, location = mu, scale = sig.sq)
    expect_gt(kt$p.value, 0.01)
  })
}

test_that('Student t draws the right number of variates', {
  ts <- rst(n = 243, df = 3)
  expect_equal(243, length(ts))
})
