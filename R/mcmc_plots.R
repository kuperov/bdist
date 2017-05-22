#' Generate MCMC trace plots based on a data.frame of draws.
#'
#' @param draws data frame of draws, one column per variable
#' @param burn.in number of burn-in draws (default 0)
#' @export
trace.plot <- function(draws, burn.in=0) {
  M <- nrow(draws) - burn.in  # used draws
  K <- ncol(draws)  # variables
  par(mfrow=c(ceiling(K/2),2))  # two columns of plot
  for (i in 1:K) {
    plot((-burn.in+1):M, draws[, i],
         type='n', main=paste0('Trace plot: ',colnames(draws)[i]),
         xlab=NA, ylab=NA)
    lines((-burn.in+1):0, draws[1:burn.in, i], col='gray')
    lines(1:M, draws[(burn.in + 1):(burn.in + M), i])
    abline(v=0)
  }
}

#' Marginal density plots, optionally discarding burn-in draws
#'
#' @param draws data frame of draws, one column per variable
#' @param burn.in number of draws to drop (default 0)
#' @export
marginal.plots <- function(draws, burn.in=0) {
  K <- ncol(draws)
  par(mfrow=c(ceiling(K/2),2))  # 2 columns of plots
  for (i in 1:K) {
    dens <- density(draws[(1+burn.in):nrow(draws), i])
    plot(dens, main=paste0('Marginal density: ',colnames(draws)[i]))
  }
}

#' Acceptance rate plot
#'
#' @param accepted vector of binary variable for acceptances (TRUE or 1=accepted,
#'                 FALSE or 0=rejected)
#' @param burn.in burn-in draws to drop
#' @export
acceptance.rate.plot <- function(accepted, burn.in=0) {
  stopifnot(length(accepted)>0, burn.in>=0)
  M <- length(accepted) - burn.in
  par(mfrow=c(1,1))
  avg.pc <- 100*mean(accepted[(1+burn.in):(M+burn.in)])
  acc.rate <- stats::filter(accepted, rep(1/13, 13))*100
  plot((1-burn.in):M, acc.rate, type='n',
       main = 'M-H Proposal Acceptance Rate', xlab='Iteration',
       ylab = 'Acceptance rate (%)')
  if (burn.in > 0) {
    lines((1-burn.in):0, acc.rate[1:burn.in], col='gray')
    abline(v=0)
  }
  lines(1:M, acc.rate[(burn.in+1):(burn.in+M)])
  mtext(sprintf('13-iteration centered moving average; post burn-in overall average %2.0f%%', avg.pc), side = 3)
}
