
#' Weibull distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the weibull distribution. This is different from the standard R
#' weibull distribution, differing in its parameterization.
#'
#' The density of a weibull-distributed variable X ~ W(alpha, beta) is
#' f(t) = gamma * alpha * t^(alpha-1) * exp(-gamma * t^alpha)
#'
#' @param n number of variates to draw
#' @param q vector of quantiles
#' @param alpha alpha parameter
#' @param gamma gamma parameter
#' @export
rweibull2 <- function(n, alpha, gamma)
  qweibull2(q = runif(n), alpha = alpha, gamma = gamma)

#' @export
#' @rdname rweibull2
dweibull2 <- function(x, alpha, gamma)
  gamma * alpha * x^(alpha - 1) * exp(-gamma*x^alpha)

#' @export
#' @rdname rweibull2
qweibull2 <- function(q, alpha, gamma)
  (-log(1 - q)/gamma)^(1/alpha)

#' @export
#' @rdname rweibull2
pweibull2 <- function(q, alpha, gamma)
  1 - exp(-gamma * q^alpha)
