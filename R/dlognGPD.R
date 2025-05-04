#' Density of the lognormal-GPD mixture
#'
#' This function evaluates the lognormal-GPD mixture density function.
#' @param x vector (nx1): points where the function is evaluated.
#' @param p real, 0<p<1: prior probability 
#' @param mu real: log-mean of the truncated lognormal distribution.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xi real: shape parameter of the generalized Pareto distribution.
#' @param beta real: scale parameter of the generalized Pareto distribution.
#' @return ydens (n x 1) vector: numerical values of the
#' lognormal - generalized Pareto mixture at x.
#' @export
#' @examples
#' ydens <- dlognGPD(seq(0,20,length.out=500),.9,0,1,0.5,2)
#'
#' @importFrom Rdpack reprompt

dlognGPD <- function(x,p,mu,sigma,xi,beta)
{
ydens <- p * dlnorm(x,mu,sigma) + (1-p) * evd::dgpd(x,0,beta,xi)
return(ydens)
}