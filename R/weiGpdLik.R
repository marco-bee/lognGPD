#' Weighted GPD log-likelihood
#'
#' This function evaluates the zero-mean generalized Pareto log-likelihood
#' function computed with weighted observations.
#' @param x numerical vector (2x1): values of the parameters \eqn{\xi}
#' and \eqn{\beta}. 
#' @param y numerical vector (nx1): observed data.
#' @param post numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @return llik real number: numerical value of the log-likelihood
#' function
#' @export
#' @examples
#' y <- rlognGPD(100,.9,0,1,0.5,2)
#' x0 <- c(.7,.2,1.3,.8,1.7)
#' res <- EMlogngpdmix(x0, y, 1000)
#' llik <- weiGpdLik(c(res$beta,res$xi),y,res$post)
#'
#' @importFrom Rdpack reprompt

weiGpdLik <- function(x,y,post)
{
  # x[1] = xi, x[2]= beta  
  llik = sum((evd::dgpd(y,0,exp(x[2]),x[1], log=TRUE))*post)
  return(llik)
}
