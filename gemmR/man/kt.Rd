\name{kt}
\alias{kt}
\title{quick Kendall's tau}
\usage{
  kt(arr1, arr2, length)
}
\arguments{
  \item{arr1}{numeric vector}

  \item{arr2}{numeric vector}

  \item{length}{length of arr1 or arr1, (both must be the same length).}
}
\value{
  \item{kendallNlogN}{Kendall's tau between \code{arr1} and \code{arr2}}
}
\description{
  This is a better-scaling implementation of Kendall's tau than is used by
  \code{\link{cor.test, method = "kendall"}}.
}
