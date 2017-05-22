# Truncated and multivariate normal distribution functions

#' Draw 1 variate from truncated normal distribution.
#' Runs in Op(1) time and space. Stack can overrun if range too far from draws.
#'
#' @param lower lower truncation bound
#' @param upper upper truncation bound
#' @param mean underlying distribution mean
#' @param sd underlying distribution sd
#' @export
rnorm1_trunc <- function(lower, upper, mean, sd) {
  stopifnot(lower < upper, sd > 0, !is.nan(mean), !is.nan(sd))
  y <- rnorm(1, mean = mean, sd = sd)
  if (y < lower || y > upper)
    rnorm1_trunc(lower, upper, mean, sd) # recursively try again
  else
    y
}

#' Density of a truncated normal distribution
#'
#' @param x quantile at which to evaluate density
#' @param lower lower truncation bound
#' @param upper upper truncation bound
#' @param mean underlying distribution mean
#' @param sd underlying distribution sd
#' @param log if true, return log density
#' @export
dnorm_trunc <- function(x, lower = -Inf, upper = Inf, mean = 0, sd = 1, log = F) {
  trunc.fac <- 1 - (pnorm(q = lower, mean = mean, sd = sd)
                    + pnorm(q = upper, mean = mean, sd = sd, lower.tail = F))
  stopifnot(lower < upper, all(trunc.fac > 0))
  if (log) {
    ifelse(x < lower | x > upper, -Inf,
           dnorm(x, mean = mean, sd = sd, log = T) - log(trunc.fac))
  }
  else {
    ifelse(x < lower | x > upper, 0,
           dnorm(x, mean = mean, sd = sd, log = F) / trunc.fac)
  }
}

#' Multivariate normal distribution
#'
#' @param x vector of quantiles
#' @param n number of variates to draw
#' @param mu mean, a d-element vector
#' @param Sigma symmetric, positive definite d*d variance-covariance matrix
#' @param tol tolerance for checking positive definiteness
#' @references Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2014). Bayesian data analysis (3E). Boca Raton, FL, USA: Chapman & Hall/CRC.
#' @export
mvdnorm <- function(x, mu, Sigma, tol=1e-6) {
  d <- length(mu)
  stopifnot(length(mu) == d, dim(Sigma) == c(d,d))
  evals <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values
  if (any(evals < -evals[1]*tol)) stop('Sigma must be positive definite')
  if (class(x) == 'matrix') {
    # must be n*d, one set of coordinates per row
    if (ncol(x) != d) stop('X must have d columns')
    n <- nrow(x)
    mu.matrix <- matrix(rep(mu, n), ncol = n)
    dens <- (2*pi)^(-d/2) /
      sqrt(prod(evals)) *
      exp(-1/2 * diag((x-t(mu.matrix)) %*% solve(Sigma) %*% t(x-t(mu.matrix))))
    drop(dens)
  }
  else {
    dens <- (2*pi)^(-d/2) /
      sqrt(prod(evals)) *
      exp(-1/2 * t(x - mu) %*% solve(Sigma) %*% (x - mu))
    drop(dens)
  }
}

#' @export
#' @rdname mvdnorm
mvrnorm <- function(n, mu, Sigma, tol=1e-6) {
  d <- length(mu)
  stopifnot(dim(Sigma) == c(d,d))
  evals <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values
  if (any(evals < -evals[1]*tol)) stop('Sigma must be positive definite')
  A <- chol(Sigma)
  mu.matrix <- matrix(rep(mu, n), n, byrow = T)
  z <- matrix(rnorm(d*n), n, byrow = T)
  mu.matrix + z %*% A
}
