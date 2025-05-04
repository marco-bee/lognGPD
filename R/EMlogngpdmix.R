#' Mixture estimation via EM
#'
#' This function estimates a static lognormal - generalized Pareto mixture
#' by means of the EM algorithm.
#' @param x0 numerical vector (5x1): initial values of the parameters p,
#' \eqn{\mu}, \eqn{\sigma}, \eqn{\xi}, \eqn{\beta}. 
#' @param y vector: observed data.
#' @param maxiter positive integer: maximum number of iterations of the EM algorithm.
#' @return a list with the following elements is returned:
#' "p" = estimated value of p,
#' "post" = posterior probabilities of all observations,
#' "mu" = estimated value of \eqn{\mu},
#' "sigma " = estimated value of \eqn{\sigma},
#' "xi" = estimated value of \eqn{\xi},
#' "beta" = estimated value of \eqn{\beta},
#' "loglik" = maximimzed log-likelihood,
#' "nit" = number of iterations
#' @export
#' @examples
#' y <- rlognGPD(100,.9,0,1,0.5,2)
#' x0 <- c(.7,.2,1.3,.8,1.7)
#' res <- EMlogngpdmix(x0, y, 1000)
#'
#' @importFrom Rdpack reprompt

EMlogngpdmix <- function(x0, y, maxiter)
{
  prior <- c(x0[1],1-x0[1])
  mu <- x0[2]
  sigma <- x0[3]
  xi <- x0[4]
  logbeta <- log(x0[5])
  N <- length(y)                                     # number of observations
  eps <- 10^(-5)                                   # convergence criterion
  change <- 1                           # initial test value for convergence
  nit <- 0                                    # initialize iteration counter
  param <- c(prior, x0[2], x0[3], x0[4], x0[5])            # arrange all parameters in a vector
  post <- matrix(0, N, 2)		# open matrix for posterior probabilities
  while (change > eps && nit <= maxiter)             # start iterations
  {
    options(digits = 6)	    # format for the display of numerical results
    parold <- param                               # store old parameter values
    f1 <- dlnorm(y,mu,sigma)
    f2 <- evd::dgpd(y,0,exp(logbeta),xi)
    f <- prior[1] * f1 + prior[2] * f2       # mixture density
    loglik <- sum(log(f))                  # evaluate log-likelihood function
    
    post[,1] <- prior[1] * f1 / f
    post[,2] <- prior[2] * f2 / f
    prior <- colMeans(post)              # M-step: prior probabilities
    
    mu <- post[,1] %*% log(y) / (N * prior[1])
    sigma <- sqrt(t(post[,1]) %*% (log(y)-as.vector(mu))^2 / (N*prior[1]))
    gpdpars = optim(c(xi,logbeta),weiGpdLik,gr=NULL,y,post[,2],control=list(fnscale=-1))
    xi = gpdpars$par[1]
    logbeta = gpdpars$par[2]

    param <- c(prior, mu, sigma, xi, logbeta)              # arrange all parameters in a vector
    change <- max(abs(param - parold))           # test value for convergence
    nit <- nit + 1                               # increase iteration counter
  }
  results <- list("p" = prior, "post" = post, "mu" = mu, "sigma " = sigma, "xi" = xi, "beta " = exp(logbeta), "loglik" = loglik, "nit" = nit)
  results
}
