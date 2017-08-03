#------ Univariate Student t -------

#' Student t distribution with location mu, scale sigma, and nu degrees of freedom.
#'
#' The pdf of the Student t distribution is given in Gelman et al (2014, p.578):
#' \deqn{t(x | \mu, \sigma, \nu) = \frac{\gamma\left(\frac{\nu+1}{2}\right)}{\gamma(\frac{\nu}{2})\sqrt{\nu\pi}\sigma}\left(1 + \frac{1}{\nu}\left(\frac{x-\mu}{\sigma}\right)^2\right)^{-\frac{\nu+1}{2}}}
#'
#' @param x vector of quantiles
#' @param q vector of quantiles
#' @param df degrees of freedom
#' @param location mu parameter
#' @param scale parameter (can be written sigma or sigma^2; this is 'sigma^2' in the above expression)
#' @param n number of random deviates to draw
#' @param log return logarithm of value if true
#' @return `dst` gives the density and `rst` generates random deviates.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2014). Bayesian data analysis (3E). Boca Raton, FL, USA: Chapman & Hall/CRC.
#' @export
dst <- function(x, df, location=0, scale=1, log=FALSE) {
  stopifnot(df > 0, scale > 0)
  #lt <- lgamma((df+1)/2) - lgamma(df/2) - lgamma(1/2) - 1/2*log(df) - log(scale) -(df+1)/2 * log(1 + 1/df*((x-location)/scale)^2)
  #if (log) lt else exp(lt)
  x.std <- (x - location)/scale
  if (log)
    dt(x.std, df = df, log = T) - log(scale)
  else
    dt(x.std, df = df) / scale
}

#' @rdname dst
#' @export
rst <- function(n, df, location=0, scale=1) {
  stopifnot(df > 0, scale > 0)
  #location + rnorm(n=n, sd=scale)/sqrt(rchisq(n=n, df=df)/df)
  location + rt(n, df = df) * scale
}

#' @rdname dst
#' @export
pst <- function(q, df, location=0, scale=1) {
  stopifnot(df > 0, scale > 0)
  pt((q-location)/scale, df = df)
}

#------ Multivariate Student t -------

#' Multivariate student t distribution
#'
#' The pdf of the Student t distribution is given in Gelman et al (2014, p.578):
#' \deqn{(x | \mu, \sigma, \nu) = \frac{\gamma\left(\frac{\nu+d}{2}\right)}{\gamma(\frac{\nu}{2})\nu^{d/2}\pi^{d/2}}\left|\Sigma\right|^{-1/2}\left(1 + \frac{1}{\nu}(x-\mu)'\Sigma^{-1}(x-\mu)\right)^{-\frac{\nu+d}{2}}}
#'
#' @param x point at which to calculate ordinate of multivariate t density
#' @param nu degrees of freedom, a scalar
#' @param mu location, a numeric vector of length d
#' @param Sigma scale, a symmetric dxd positive definite matrix
#' @param n number of random deviates to draw
#' @param tol tolerance for checking positive definiteness
#' @return `mvdst` gives the density and `mvrst` generates random deviates.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2014). Bayesian data analysis (3E). Boca Raton, FL, USA: Chapman & Hall/CRC.
#' @name mvst
NULL

#' @rdname mvst
#' @export
dmvst <- function(x, nu, mu, Sigma, tol=1e-6) {
  d <- length(x)
  evals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(length(nu) == 1, length(mu) == d, dim(Sigma) == c(d,d))
  if (any(evals < -evals[1]*tol)) stop('Sigma must be positive definite')
  dens <- gamma((nu+d)/2) *
    prod(evals)^(-1/2) *
    (1+t(x-mu)%*%solve(Sigma)%*%(x-mu)/nu)^(-(nu+d)/2) /
    (gamma(nu/2)*nu^(d/2)*pi^(d/2))
  drop(dens)
}

#' @rdname mvst
#' @export
rmvst <- function(n, nu, mu, Sigma, tol=1e-6) {
  d <- length(mu)
  evals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(length(nu) == 1, dim(Sigma) == c(d,d))
  if (any(evals < -evals[1]*tol)) stop('Sigma must be positive definite')
  z <- mvrnorm(n=n, mu=rep(0,d), Sigma=Sigma, tol=tol)
  x <- rchisq(n=n, df=nu)
  mu.matrix <- matrix(rep(mu, n), byrow = T, ncol = d)
  mu.matrix + z * sqrt(nu/x)
}
