# generalized double pareto distribution

#' Generalized double Pareto distribution.
#'
#' If X ~ gdP(alpha, xi), then it has density
#'
#' f(X|alpha,xi) = 1/(2*xi) * (1 + abs(X)/(alpha*xi))^(-alpha-1)
#'
#' @param xi scale parameter, where xi > 0
#' @param alpha shape parameter, where alpha > 0
#' @param n number of variates to simulate
#' @param x vector of quantiles
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @references Armagan, Artin, David B. Dunson, and Jaeyong Lee. "Generalized
#'             double Pareto shrinkage." Statistica Sinica 23.1 (2013): 119.
#' @export
#' @rdname gdp
dgdp <- function(x, xi, alpha) {
  stopifnot(xi > 0, alpha > 0)
  1/(2*xi) * (1 + abs(x)/(alpha*xi))^(-alpha-1)
}

#' @rdname gdp
#' @export
rgdp <- function(n, xi, alpha) {
  stopifnot(xi > 0, alpha > 0)
  # the gdp is a three level scale mixture of normals
  eta <- xi*alpha
  lambda <- rgamma(n, shape = alpha, rate = eta)
  tau <- rexp(n, lambda^2/2)
  rnorm(n, 0, sd=sqrt(tau))
}

#' @rdname gdp
#' @export
pgdp <- function(q, xi, alpha) {
  stopifnot(xi > 0, alpha > 0)
  ifelse(q >= 0,
         1/2 * (2 - ((alpha * xi)/(q + alpha * xi))^alpha),
         1/2 * ((alpha * xi)/(-q + alpha * xi))^alpha)
}

#' @rdname gdp
#' @export
qgdp <- function(p, xi, alpha) {
  stopifnot(xi > 0, alpha > 0)
  ifelse(p >= 1/2,
         alpha * xi * ( (2-2*p)^(-1/alpha) - 1 ),
         alpha * xi * ( 1 - (2*p)^(-1/alpha) ) )
}
