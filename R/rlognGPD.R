#' Simulation of the lognormal-GPD mixture
#'
#' This function simulates the lognormal-GPD mixture.
#' @param n positive integer: number of observations sampled.
#' @param p real, 0<p<1: prior probability 
#' @param mu real: log-mean of the truncated lognormal distribution.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xi real: shape parameter of the generalized Pareto distribution.
#' @param beta real: scale parameter of the generalized Pareto distribution.
#' @return ysim (n x 1) vector: n random numbers from the
#' lognormal - generalized Pareto mixture.
#' @export
#' @examples
#' ysim <- rlognGPD(100,.9,0,1,0.5,2)
#'
#' @importFrom Rdpack reprompt

rlognGPD <- function(n,p,mu,sigma,xi,beta)
{
  z <- rbinom(n,1,p)
  n1 <- sum(z)
  n2 <- n - n1
  y1 <- rlnorm(n1,mu,sigma)
  y2 <- evd::rgpd(n2,0,beta,xi)
  ysim <- c(y1,y2)
  return(ysim)
}
