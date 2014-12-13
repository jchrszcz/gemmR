\name{gemm}
\alias{gemm}
\alias{gemm.default}
\alias{gemm.formula}
\alias{print.gemm}
\alias{plot.gemm}
\alias{summary.gemm}
\alias{print.summary.gemm}
\alias{nobs.gemm}
\alias{deviance.gemm}
\alias{logLik.gemm}
\title{Fit General Monotone Models.}
\usage{
\method{gemm}{formula}(formula, data=list(), ...)

\method{gemm}{default}(x, k.pen, parallel = FALSE, n.chains = 4, fit.metric = "bic", ...)
}
\arguments{
  \item{formula}{an object of the class \code{\link{formula}}.}
  
  \item{data}{optional data frame, list, or environment.}
  
  \item{x}{must be data frame, first column is treated as dependent variable.
  Not needed if formula is supplied.}
  
  \item{k.pen}{vector of integers giving penalty equivalent to main effects for
  any interaction terms. Calculated by \code{\link{gemm.formula}} or defaults to
  no special interaction weighting.}
  
  \item{parallel}{boolean flag to run chains in parallel over available cores.}
  
  \item{n.chains}{number or replications of the sampling process, used to
  assess various starting conditions}

  \item{fit.metric}{Value to optimize with genetic algorithm. Currently accepts ``bic'', ``aic'', and ``tau''.}

  \item{...}{Additional arguments to be passed to lower level fitting functions (see below).}
}
\description{
  \code{gemm} is used to fit general monotone models. By default, the function
  will generate metric weights that minimize rank order inversion between the
  model predictions and a response variable, subject to a parsimony correction.
  Optional argument passed to \code{\link{gemmEst}}.
}
\details{
  Models for \code{gemm} are constructed with syntax similar to
  \code{\link{lm}}. By default, \code{gemm} will use random search to minimize
  penalized rank order inversions between model predictions and a response
  variable. This is accomplished by generating candidate weights using
  \code{\link{genAlg}} and deriving a BIC based on transformed Kendall's
  \code{tau} using \code{\link{gemmFitRcppI}}. \code{gemm} may take some time
  to run depending on the complexity of the model. Methods for standard
  goodness-of-fit functions are also available and are run by
  \code{\link{summary.gemm}}.
}
\value{
  A list with class "\code{gemm}" containing the following components:
  
  \item{date}{system time and date for model completion.}
  
  \item{call}{the matched call.}
  
  \item{coefficients}{matrix of best weights with one row for each chain.}
  
  \item{fitted.values}{model predictions for criterion generated from weights
  associated with best chain.}
  
  \item{residuals}{metric values for response minus fitted.}
  
  \item{rank.residuals}{rank response minus rank criterion.}
  
  \item{bic}{vector of Bayesian Information Criteria for estimation sample
  of each chain.}
  
  \item{aic}{vector of Aikaike Information Criteria for estimation sample
  of each chain.}
  
  \item{r}{vector of Pearson's \code{r} values for estimation sample of
  each chain.}
  
  \item{tau}{vector of Kendall's \code{tau} for estimation sample of each
  chain.}
  
  \item{metric.betas}{regression weights derived using \code{\link{lm}}.}
  
  \item{p.vals}{p-values associated with ordinary least squares regression
  weights.}
  
  \item{model}{data frame including modeled data.}
  
  \item{fit.metric}{sorting metric used.}
  
  \item{cross.val.bic}{vector of Bayesian Information Criteria for
  cross-validation sample of each chain.}
  
  \item{cross.val.aic}{vector of Aikaike Information Criteria for
  cross-validation sample of each chain.}

  \item{cross.val.r}{vector of Pearson's \code{r} values for cross-validation
  sample of each chain.}
  
  \item{cross.val.tau}{vector of Kendall's \code{tau} for cross-validation
  sample of each chain.}
  
  \item{converge.fit.metric}{matrix of "\code{fit.metric}" values with
  generations by chains.}
  
  \item{converge.beta}{matrix derived weights for each generation within each
  chain, column for each predictor.}
  
  \item{converge.r}{generations by chains matrix of Pearson's \code{r}.}
  
  \item{formula}{\code{\link{formula}} object used.}
}
\seealso{
  \code{\link{genAlg}} for search, \code{\link{gemmFitRcppI}} for fitting
  routine, \code{\link{tauTest}} for \code{O*(N log N)} scale Kendall's \code{tau},
  \code{\link{convergencePlot}} for the optional plot pane when
  \code{check.convergence = TRUE}. \code{\link{gemmEst}} may be useful in some
  circumstances when formula input and interaction terms are not needed.
}
\examples{
  \dontrun{
    data(culture)
    gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars, check.convergence = TRUE)
    print(gemm.model)
    plot(gemm.model)
  }
}
\references{
  Dougherty, M. R., & Thomas, R. P. (2012). Robust decision making in a
  nonlinear world. \emph{Psychological review, 119}(2), 321.
}
\keyword{
  ordinal, regression
}
