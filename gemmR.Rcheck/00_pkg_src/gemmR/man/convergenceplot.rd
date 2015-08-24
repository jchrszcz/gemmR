\name{convergencePlot}
\alias{convergencePlot}
\title{convergencePlot function for visualizing genetic algorithm in gemmR}
\usage{
  convergencePlot(beta, fit.metric, ...)
}
\arguments{
  \item{beta}{matrix of BIC values for the best model in
  each generation.}
  
  \item{fit.metric}{selected metric for a given model, specified when
  \code{\link{gemm}} is run.}

  \item{...}{Additional arguments to be passed to lower level plotting functions.}
}
\description{
  This function plots convergence for gemm objects that
  were called with the "check.convergence" argument.
}
\details{
  Generally does not work well with \code{n.chains} < 10.
}
\examples{
  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + hp, data = mtcars, check.convergence = TRUE,
    seed.metric = FALSE, n.chains = 3, n.gens = 3, n.beta = 200)
  with(gemm.model, convergencePlot(converge.fit.metric, fit.metric))
}