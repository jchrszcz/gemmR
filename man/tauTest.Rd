\name{tauTest}
\alias{tauTest}
\title{quick Kendall's tau}
\usage{
  tauTest(arr1, arr2, length)
}
\arguments{
  \item{arr1}{numeric vector}

  \item{arr2}{numeric vector}

  \item{length}{length of arr1 or arr1, (both must be the same length).}
}
\value{
  Returns a list with the following components:

  \item{tau-a}{Kendall's tau-a, which includes tied pairs of all types in the denominator.}

  \item{tau-b}{Kendall's tau-b, which does not penalize for ties.}

  \item{n.pairs}{Number of paired comparisons.}

  \item{n.ties.1}{Number of ties in arr1.}

  \item{n.ties.2}{Number of ties in arr2.}

  \item{n.ties.both}{Number of ties for both arr1 and arr2.}

  \item{n.dis}{Total disconcordances for all paired comparisons.}

  \item{n.con}{Total concordances for all paired comparisons.}
}
\description{
  This is a better-scaling implementation of Kendall's tau than is used by
  \code{\link{cor.test}}.
}
