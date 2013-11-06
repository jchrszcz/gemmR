\name{convergencePlot}
\alias{convergencePlot}
\title{convergencePlot function for visualizing genetic algorithm in gemmR}
\usage{
  convergencePlot(beta, ...)
}
\arguments{
  \item{beta}{Matrix of BIC values for the best model in
  each generation.}
}
\description{
  This function plots convergence for gemm objects that
  were called with the "check.convergence" argument.
  Generally doesn't work well with n.data.gen < 10.
}

