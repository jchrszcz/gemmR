\name{kCorFact}
\alias{kCorFact}
\title{Compute BIC penalty for factors and interaction terms}
\usage{
  kCorFact(k.pen, beta.vecs)
}
\arguments{
  \item{k.pen}{matrix identifying which terms are relevant to a given input.}

  \item{beta.vecs}{matrix of candidate beta weights.}
}
\value{
  \item{out}{vector of effective terms included by virture of a non-zero weight, accounting for interaction and factor grouings.}
}
\description{
  This is used by \code{\link{gemmEst}} to calculate the appropriate penalty for BIC and AIC calculations.
  This algorithm generates the \code{k} term in the calculation of information criteria.
}