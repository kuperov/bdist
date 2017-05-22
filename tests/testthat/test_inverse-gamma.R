context('gamma-ish distributions')

set.seed(1234) # can comment out for non-automated tests

params <- data.frame(shape=c(1,6,2,3,32,2,40), df=c(2,7,4,2,1,25,35))

integral.tol <- 1e-4

for (i in 1:nrow(params)) {
  shape <- params$shape[i]
  df <- params$df[i]

  test_that(sprintf('Inverse gamma density integrates to 1 for shape = %f, scale = %f', shape, df), {
    igral <- integrate(function(x) dinvgamma(x, shape = shape, scale = df), lower=0, upper=Inf)
    expect_equal(1, igral$value)
  })

  test_that(sprintf('Inverse root gamma density integrates to 1 for shape = %f, df = %f', shape, df), {
    igral <- integrate(function(x) dinvrootgamma(x, shape = shape, df = df), lower=0, upper=Inf)
    expect_equal(1, igral$value)
  })

  test_that(sprintf('Inverse root gamma density integrates to cdf for shape = %f, df = %f', shape, df), {
    limits <- seq(1,5)
    lapply(limits, function(lim) {
      igral <- integrate(function(x) dinvrootgamma(x, shape = shape, df = df), lower=0, upper=lim)
      p <- pinvrootgamma(lim, shape = shape, df = df)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('Inverse gamma density integrates to cdf for shape = %f, scale = %f', shape, df), {
    limits <- seq(1,5)
    lapply(limits, function(lim) {
      igral <- integrate(function(x) dinvgamma(x, shape = shape, scale = df), lower=0, upper=lim)
      p <- pinvgamma(lim, shape = shape, scale = df)
      expect_lt(abs(igral$value - p), integral.tol)
    })
  })

  test_that(sprintf('KS test for inverse root gamma draws & cdfs (shape = %.1f, df = %.1f)', shape, df), {
    draws <- rinvrootgamma(1e4, shape = shape, df = df)
    kt <- ks.test(draws, 'pinvrootgamma', shape = shape, df = df)
    expect_gt(kt$p.value, 0.01)
  })

  test_that(sprintf('KS test for inverse gamma draws & cdfs (shape = %.1f, rate = %.1f)', shape, df), {
    draws <- rinvgamma(1e4, shape = shape, scale = df)
    kt <- ks.test(draws, 'pinvgamma', shape = shape, scale = df)
    expect_gt(kt$p.value, 0.01)
  })
}
