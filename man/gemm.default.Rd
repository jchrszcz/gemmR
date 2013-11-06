\name{gemm.default}
\alias{gemm.default}
\title{Default method for gemm}
\usage{
  gemm.default(x, k.pen = k.pen, ...)
}
\arguments{
  \item{x}{dataframe to be passed to \code{\link{gemmEst}}.
  First column is predicted variable.}

  \item{k.pen}{Penalty term, necessary for appropriate
  penalty related to interaction terms.}

  \item{...}{Other arguments passed to
  \code{\link{gemmEst}}.}
}
\description{
  This function is the default method for gemm. It calls
  \code{\link{gemmEst}}, assigns a "gemm" class to the
  resultant list object, and returns that list.
}

