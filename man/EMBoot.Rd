% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMBoot.R
\name{EMBoot}
\alias{EMBoot}
\title{Bootstrap standard errors for the MLEs of a lognormal-GPD mixture}
\usage{
EMBoot(x, x0, y, maxiter)
}
\arguments{
\item{x}{list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
of a single dataset.}

\item{y}{numerical vector: observed sample.}

\item{maxiter}{non-negative integer: maximum number of iterations of the EM algorithm.}
}
\value{
Estimated parameters obtained from a bootstrap sample.
}
\description{
This function draws a bootstrap sample and uses it to estimate the parameters of a lognormal-Pareto mixture distribution. Since this is typically called by LPfitEM, see the help of LPfitEM for examples.
}
\details{
At each bootstrap replication, the mixture is estimated via the EM algorithm.
}
