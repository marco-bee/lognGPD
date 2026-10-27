#' Cumulative distribution function of the lognormal-GPD mixture
#'
#' This function evaluates the lognormal-GPD mixture cumulative distribution function.
#' @param x vector (nx1): points where the function is evaluated.
#' @param p real, 0<p<1: prior probability 
#' @param mu real: log-mean of the truncated lognormal distribution.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xi real: shape parameter of the generalized Pareto distribution.
#' @param beta positive real: scale parameter of the generalized Pareto distribution.
#' @return ycdf (n x 1) vector: numerical values of the
#' lognormal - generalized Pareto cdf at x.
#' @import stats evd
#' @export
#' @examples
#' ycdf <- plognGPD(seq(0,20,length.out=500),.9,0,1,0.5,2)
#'
#' @importFrom Rdpack reprompt

plognGPD <- function(x,p,mu,sigma,xi,beta)
{
  ycdf <- p * plnorm(x,mu,sigma) + (1-p) * evd::pgpd(x,0,beta,xi)
  return(ycdf)
}
