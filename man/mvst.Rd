% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/student.R
\name{mvst}
\alias{mvst}
\alias{dmvst}
\alias{rmvst}
\title{Multivariate student t distribution}
\usage{
dmvst(x, nu, mu, Sigma, tol = 1e-06)

rmvst(n, nu, mu, Sigma, tol = 1e-06)
}
\arguments{
\item{x}{point at which to calculate ordinate of multivariate t density}

\item{nu}{degrees of freedom, a scalar}

\item{mu}{location, a numeric vector of length d}

\item{Sigma}{scale, a symmetric dxd positive definite matrix}

\item{tol}{tolerance for checking positive definiteness}

\item{n}{number of random deviates to draw}
}
\value{
`mvdst` gives the density and `mvrst` generates random deviates.
}
\description{
The pdf of the Student t distribution is given in Gelman et al (2014, p.578):
\deqn{(x | \mu, \sigma, \nu) = \frac{\gamma\left(\frac{\nu+d}{2}\right)}{\gamma(\frac{\nu}{2})\nu^{d/2}\pi^{d/2}}\left|\Sigma\right|^{-1/2}\left(1 + \frac{1}{\nu}(x-\mu)'\Sigma^{-1}(x-\mu)\right)^{-\frac{\nu+d}{2}}}
}
\references{
Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2014). Bayesian data analysis (3E). Boca Raton, FL, USA: Chapman & Hall/CRC.
}
