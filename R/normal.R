# Truncated and multivariate normal distribution functions

#' Truncated normal distribution.
#'
#' Density, distribution function, quantile function, and random generation
#' for the normal distribution truncated at (lower, upper).
#'
#' @param x vector of quantiles
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of variates to generate
#' @param mean underlying distribution mean
#' @param sd underlying distribution sd
#' @param lower lower truncation bound
#' @param upper upper truncation bound
#' @param log return logarithm
#' @export
#' @rdname norm_trunc
rnorm_trunc <- function(n, mean, sd, lower, upper) {
  stopifnot(n > 0, lower < upper, sd > 0, !is.nan(mean), !is.nan(sd))
  qnorm_trunc(runif(n), mean = mean, sd = sd, lower = lower, upper = upper)
}

#' @export
#' @rdname norm_trunc
dnorm_trunc <- function(x, lower = -Inf, upper = Inf, mean = 0, sd = 1, log = F) {
  trunc.fac <- pnorm(upper, mean = mean, sd = sd) -
    pnorm(lower, mean = mean, sd = sd)
  stopifnot(lower < upper, all(trunc.fac > 0))
  if (log) {
    ifelse(x < lower | x > upper,
           -Inf,
           dnorm(x, mean = mean, sd = sd, log = T) - log(trunc.fac))
  }
  else {
    ifelse(x < lower | x > upper,
           0,
           dnorm(x, mean = mean, sd = sd, log = F) / trunc.fac)
  }
}

#' @export
#' @rdname norm_trunc
pnorm_trunc <- function(q, mean, sd, lower, upper) {
  ifelse(q < lower, 0,
         ifelse(q > upper, 1,
                (pnorm(q, mean, sd) - pnorm(lower, mean, sd))/
                  (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
                ))
}

#' @export
#' @rdname norm_trunc
qnorm_trunc <- function(p, mean, sd, lower, upper) {
  ifelse(p < 0 | p > 1, NA,
  ifelse(p == 0, -Inf,
  ifelse(p == 1, upper,
         qnorm(p * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) +
               pnorm(lower, mean, sd),
               mean, sd)
         )))
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
