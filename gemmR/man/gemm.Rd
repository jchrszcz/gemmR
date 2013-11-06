\name{gemm}
\alias{gemm}
\title{gemm}
\usage{
  gemm(x, ...)
}
\arguments{
  \item{x}{Formula input for gemm procedure.}

  \item{...}{Other arguments to be passed to
  \code{\link{gemm.formula}} and \code{\link{gemmEst}}.}
}
\description{
  This function is called by \code{\link{gemm.default}} and
  runs a handful of subordinate functions to produce a gemm
  object. Takes data and, over successive replications,
  uses geneticAlgorithm to generate candidate beta vectors,
  calculates ordinal model fit using these betas, and
  produces an output that reports weights and fit
  statistics for best models at each generation,
  (optionally, for cross-validation as well).
}

