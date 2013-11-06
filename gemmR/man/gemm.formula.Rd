\name{gemm.formula}
\alias{gemm.formula}
\title{Function method for gemm}
\usage{
  gemm.formula(formula, data = list(), ...)
}
\arguments{
  \item{formula}{Linear model equivalent to be fit by
  gemm.}

  \item{data}{Dataframe referred to by formula.}

  \item{...}{Other arguments to be passed to
  \code{\link{gemmEst}}.}
}
\description{
  This is the function method for gemm. It calls
  \code{\link{gemm.default}} after creating the design
  matrix and calculating the appropriate penalty term
}

