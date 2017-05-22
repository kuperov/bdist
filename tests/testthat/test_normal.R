# test multivariate normal distribution functions

library(cubature)

tol <- 1e-4
test_that('2D multivariate normal integrates to unity', {
  # we'll use 50 as 'infinity'
  mu <- c(1,2)
  Sigma <- matrix(c(3,0.5,0.5,6), 2)
  f <- function(x)
    matrix(mvdnorm(t(x), mu=mu, Sigma=Sigma), ncol=1)
  volume <- pcubature(f, lowerLimit=c(-50, -50), upperLimit=c(50, 50), vectorInterface = TRUE)
  expect_true(abs(volume$integral - 1) < tol)
})

