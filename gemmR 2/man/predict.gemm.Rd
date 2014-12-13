\name{predict.gemm}
\alias{predict.gemm}
\title{Predict method for general monotone models.}
\usage{
\method{predict}{gemm}(object, newdata = NULL, \dots, tie.struct = FALSE)
}
\arguments{
  \item{object}{an object of the class "\code{\link{gemm}}".}
  
  \item{newdata}{optional list of new data.}

  \item{\dots}{further arguments passed to or from other methods.}
  
  \item{tie.struct}{option to return corrects, incorrects, and ties based on
  fitted values.}
}
\description{
  \code{predict.gemm} is used to get fitted values from "\code{gemm}" objects.
  Can also be used to get correct, incorrect, and ties for all paired
  comparisons of cases in model.
}
\details{
  Will eventually include functionality to generate estimates based on new data. 
}
\value{
  Output depends on \code{tie.struct} argument. Output is either:
  
  \item{y}{matrix of fitted values or a data.frame of correct, incorrect, and tied pairs..}
}
\seealso{
  \code{\link{gemm}}
}
\examples{
  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars, n.chains = 3, n.gens = 3, n.beta = 200)
  predict(gemm.model, tie.struct = TRUE)
}
\keyword{
  ordinal, regression
}
