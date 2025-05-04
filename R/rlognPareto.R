#' Simulation of the lognormal-Pareto spliced distribution
#'
#' This function fits the continuous and differentiable
#' version of the truncated lognormal-Pareto spliced distribution
#' proposed by Scollnik (2007).
#' @param n positive integer: number of observations sampled.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xmin positive real: scale parameter of the Pareto distribution.
#' @param alpha positive real: shape parameterof the Pareto distribution.
#' @return ysim (nreps x 1) vector: nreps random numbers from the
#' truncated lognormal-Pareto spliced distribution.
#' @details See Scollnik (2007) for details.
#' @keywords dynamic mixture; simulation.
#' @export
#' @examples
#' ysim <- rlognPareto(100,1,5,2)
#' @references{
#'   \insertRef{scoll07}{lognGPD}
#' }
#'
#'
#' @importFrom Rdpack reprompt

rlognPareto <- function(n,sigma,xmin,alphapar)
{
mu <- log(xmin)-alphapar*sigma^2      
p_num <- sqrt(2*pi)*alphapar*sigma*pnorm(alphapar*sigma)*
  exp(.5*(alphapar*sigma)^2)
p_den <- sqrt(2*pi)*alphapar*sigma*pnorm(alphapar*sigma)*
  exp(.5*(alphapar*sigma)^2)+1
p <- p_num/p_den
y1 <- EnvStats::rlnormTrunc(floor(n*p),mu,sigma,0,xmin)
y2 <- LNPar::rpareto(ceiling(n*(1-p)),xmin,alphapar)
y = c(y1,y2)
return(y)
}