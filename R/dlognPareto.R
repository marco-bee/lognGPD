#' Density of the lognormal-Pareto spliced distribution
#'
#' This function evaluates the density of the continuous and
#' differentiable version of the truncated lognormal-Pareto spliced
#' distribution proposed by Scollnik (2007).
#' @param x vector (nx1): points where the function is evaluated.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xmin positive real: scale parameter of the Pareto distribution.
#' @param alpha positive real: shape parameterof the Pareto distribution.
#' @return ysim (n x 1) vector: numerical values of the
#' truncated lognormal-Pareto spliced distribution at x.
#' @details See Scollnik (2007) for details.
#' @export
#' @examples
#' ysim <- dlognPareto(seq(0,20,length.out=500),1,5,2)
#' @references{
#'   \insertRef{scoll07}{lognGPD}
#' }
#'
#'
#' @importFrom Rdpack reprompt

dlognPareto <- function(x,sigma,xmin,alpha)
{
f <- rep(0,length(x))
r_num <- sqrt(2*pi)*alpha*sigma*pnorm(alpha*sigma)*exp(.5*alpha^2*sigma^2)
r_den <- sqrt(2*pi)*alpha*sigma*pnorm(alpha*sigma)*exp(.5*alpha^2*sigma^2)+1
r = r_num/r_den
mu = log(xmin) - alpha*sigma^2
indici1 = which(x<xmin)
x1 <- x[indici1]
f[indici1] <- r * dlnorm(x1,mu,sigma)/pnorm((log(xmin)-mu)/sigma)
indici2 <- which(x>xmin)
x2 <- x[indici2]
f[indici2] <- (1-r) * LNPar::dpareto(x2,xmin,alpha)
return(f)
}