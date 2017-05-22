# Inverse gamma distribution functions

#' Inverse gamma distribution (sometimes called 'reciprocal gamma')
#'
#' We follow Gelman et al's (2014) parameterization. If \eqn{X\sim Gamma(a,b)}, then \eqn{1/X\sim InvGamma(a,b)}.
#' Note that some authors use alternative parameterizations; see especially \link{dinvrootgamma}.
#'
#' The density function with shape \eqn{\alpha} and rate \eqn{\beta} is
#' \deqn{\frac{\beta ^{\alpha } x^{-\alpha -1} \exp\{-\frac{\beta }{x}\}}{\Gamma (\alpha )}.}
#'
#' The cdf with shape \eqn{\alpha} and scale \eqn{\beta} is
#' \deqn{\frac{\Gamma \left(\alpha ,\frac{\beta }{x}\right)}{\Gamma (\alpha )}.}
#'
#' @param x vector of quantiles
#' @param q vector of quantiles
#' @param n number of random deviates to draw
#' @param shape shape parameter alpha; must be positive
#' @param scale scale parameter beta; must be positive
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return `dinvgamma` gives the density and `rinvgamma` generates random deviates.
#' @references Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2014). Bayesian data analysis (3E). Boca Raton, FL, USA: Chapman & Hall/CRC.
#' @family gamma
#' @export
dinvgamma <- function(x, shape, scale, log=FALSE) {
  stopifnot(all(shape > 0), all(scale > 0))
  lig <- ifelse(x > 0, shape*log(scale) - lgamma(shape) - (shape+1)*log(x) - scale/x, -Inf)
  if (log) lig else exp(lig)
}

#' @rdname dinvgamma
#' @export
rinvgamma <- function(n, shape, scale) {
  stopifnot(all(shape > 0), all(scale > 0))
  1/(rgamma(n, shape=shape, rate=scale))
}

#' @rdname dinvgamma
#' @export
pinvgamma <- function(q, shape, scale, log = FALSE) {
  stopifnot(all(shape > 0), all(scale > 0))
  # see ?pgamma for the partial gamma function; note gamma(shape) cancels
  lp <- ifelse(q > 0, pgamma(scale/q, shape, lower.tail = F, log.p = T), -Inf)
  if (log) lp else exp(lp)
}

#' Inverse root gamma distribution (sometimes called 'inverse gamma')
#'
#' If X ~ InvRootGamma(scale=sigma.sq, df=nu), then
#' 1/(X^2) ~ Gamma(shape=nu/2, rate=nu*sigma.sq/2).
#'
#' The density function of the inverse-root-gamma distribution shape \eqn{\sigma^2}
#' and \eqn{\nu} degrees of freedom is
#' \deqn{p(\theta) = I(\theta>0)\frac{2}{\Gamma\left(\frac{\nu}{2}\right)}
#'       \left( \frac{\nu\sigma^2}{2}\right)^{\nu/2}\frac{1}{\theta^{\nu+1}}
#'       \exp\left\{-\frac{\nu \sigma^2}{2\theta^2}\right\}.}
#'
#' The cdf is
#' \deqn{P(x) = \frac{\Gamma \left(\frac{\nu }{2},\frac{\nu  \sigma ^2}{2 x^2}\right)}{\Gamma \left(\frac{\nu }{2}\right)}.}
#'
#' Note that some authors use alternative parameterizations; see especially \link{dinvgamma}.
#'
#' @param x vector of quantiles
#' @param q vector of quantiles
#' @param shape parameter, where shape>0
#' @param df degrees of freedom parameter, where df>0
#' @param n number of random deviates to draw
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return `dinvgamma` gives the density and `rinvgamma` generates random deviates.
#' @family gamma
#' @export
dinvrootgamma <- function(x, shape, df, log=FALSE) {
  stopifnot(all(shape > 0), all(df > 0))
  lirg <- ifelse(x>0, log(2) - lgamma(df/2) + df/2*log(df*shape/2) - (df+1)*log(x) - df*shape/(2*x^2), -Inf)
  if (log) lirg else exp(lirg)
}

#' @rdname dinvrootgamma
#' @export
rinvrootgamma <- function(n, shape, df) {
  stopifnot(all(shape > 0), all(df > 0))
  1/sqrt(rgamma(n, shape = df/2, rate = shape*df/2))
}

#' @rdname dinvrootgamma
#' @export
pinvrootgamma <- function(q, shape, df, log=FALSE) {
  stopifnot(all(shape > 0), all(df > 0))
  # see ?pgamma for the partial gamma function
  lq <- ifelse(q > 0,
               pgamma(q = df*shape/(2*q^2), shape = df/2, lower.tail=F, log.p=T),
               -Inf)
  if (log) lq else exp(lq)
}
