\name{gemmEst}
\alias{gemmEst}
\title{Fit General Monotone Models.}
\usage{
  gemmEst(input.data, output = "gemmr", n.beta = 8000, n.chains = n.chains,
  n.gens = 10, save.results = FALSE, k.pen = k.pen, seed.metric = TRUE,
  check.convergence = FALSE, roe = FALSE, fit.metric = fit.metric,
  correction = "knp", oclo = TRUE, isTauB = FALSE)
}
\arguments{
  \item{input.data}{must be data frame, first column is
  treated as dependent variable.}

  \item{output}{string argument for use in naming file
  output. \code{\link{gemmEst}} may write a .RData file in the
  current working directory each time the function is
  called.}

  \item{n.beta}{number of beta vectors to generate per
  generation.}

  \item{n.chains}{number of times the fitting process
  will be repeated.}

  \item{n.gens}{number of generations per chain.}

  \item{save.results}{logical value to determine whether
  the resulting gemm object is saved to a .RData file.}

  \item{k.pen}{penalty term for BIC, as calculated by
  \code{\link{gemm.formula}}.}

  \item{seed.metric}{logical value to control whether
  \code{\link{genAlg}} is seeded with OLS regression
  weights or random values.}

  \item{check.convergence}{logical value to indicate
  whether BIC for each generation is retained, mostly
  useful for checking performance of \code{\link{genAlg}}.}

  \item{roe}{logical value to determine whether region of
  equivalence data are retained.}
  
  \item{fit.metric}{value used to order models.}
  
  \item{correction}{placeholder for correction transformations on the fit
  statistic, (r or tau).}

  \item{oclo}{logical for \emph{Ordered Constrained Linear Optimization}. If 
  \code{TRUE} secondary sort based on \code{r} performed on matrix of beta 
  weights in order to maximize linear fit, given the ordered fits determined by 
  \code{fit.metric}.}

  \item{isTauB}{logical for whether to include ties in the denominator of the
  tau calculation.}

}
\description{
  \code{gemmEst} is called by \code{\link{gemm}} to fit general monotone models.
}
\details{
  Formula syntax and interaction penalty terms can be avoided by fitting data
  directly using \code{gemmEst}.
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
  
  \item{r}{vector of Pearson's \code{r} values for estimation sample of
  each chain.}
  
  \item{tau}{vector of Kendall's \code{tau} for estimation sample of each
  chain.}

  \item{tau.par}{vector containing the correction, incorrect, criterion
  ties, predictor ties, and both between the weighted cues and outcome used
  for model fit.}
  
  \item{metric.betas}{regression weights derived using \code{\link{lm}}.}
  
  \item{p.vals}{p-values associated with ordinary least squares regression
  weights.}
  
  \item{model}{data frame including modeled data.}
  
  \item{fit.metric}{sorting metric used.}
  
  \item{cross.val.bic}{vector of Bayesian Information Criteria for
  cross-validation sample of each chain.}
  
  \item{cross.val.r}{vector of Pearson's \code{r} values for cross-validation
  sample of each chain.}
  
  \item{cross.val.tau}{vector of Kendall's \code{tau} for cross-validation
  sample of each chain.}

  \item{cross.val.tau.par}{vector containing the correction, incorrect,
  criterion ties, predictor ties, and both between the weighted cues and
  outcome used for crossvalidation.}
  
  \item{converge.fit.metric}{matrix of "\code{fit.metric}" values with
  generations by chains.}
  
  \item{converge.beta}{matrix derived weights for each generation within each
  chain, column for each predictor.}
  
  \item{converge.r}{generations by chains matrix of Pearson's \code{r}.}
  
  \item{formula}{\code{\link{formula}} object used.}
}
\seealso{
  \code{\link{gemm}} for normal use. \code{\link{genAlg}} for search,
  \code{\link{gemmFitRcppI}} for fitting routine, \code{\link{tauTest}} for
  \code{O*(N log N)} scale Kendall's \code{tau}, \code{\link{convergencePlot}}
  for the optional plot pane when \code{check.convergence = TRUE}.
}
\examples{
  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars,
    check.convergence = TRUE, n.chains = 3, n.gens = 3, n.beta = 200)
  print(gemm.model)
  plot(gemm.model)
}
\references{
  Dougherty, M. R., & Thomas, R. P. (2012). Robust decision making in a
  nonlinear world. \emph{Psychological review, 119}(2), 321.
}
\keyword{
  ordinal, regression
}
