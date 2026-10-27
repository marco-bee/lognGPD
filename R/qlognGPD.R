#' Quantile function of the lognormal-GPD mixture
#'
#' This function evaluates the lognormal-GPD mixture quantile function.
#' @param w vector (nx1): points where the function is evaluated.
#' @param p real, 0<p<1: prior probability 
#' @param mu real: log-mean of the truncated lognormal distribution.
#' @param sigma positive real: log-standard deviation of the truncated 
#' lognormal distribution.
#' @param xi real: shape parameter of the generalized Pareto distribution.
#' @param beta positive real: scale parameter of the generalized Pareto distribution.
#' @return yq (n x 1) vector: numerical values of the
#' lognormal - generalized Pareto quantile function at w.
#' @import stats evd
#' @export
#' @examples
#' ycdf <- qlognGPD(seq(0,1,length.out=500),.9,0,1,0.5,2)
#'
#' @importFrom Rdpack reprompt

qlognGPD <- function(w,p,mu,sigma,xi,betap)
{
  nw <- length(w)
  yq <- rep(NA,nw)
  maxq <- max(c(qlnorm(w,mu,sigma),evd::qgpd(w,0,betap,xi))) 
  for (i in 1:nw)
  {
    yq[i] <- uniroot(function(x) plognGPD(x,p,mu,sigma,xi,betap)-w[i],
                     interval = c(0,maxq),
                     tol = 10^{-16})$root  
  }
  return(yq)
}
