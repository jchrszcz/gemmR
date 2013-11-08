\name{gemmEst}
\alias{gemmEst}
\title{gemmR controller function.}
\usage{
  gemmEst(input.data, output = "gemmr", n.beta = 8000, p.est = 1,
                    n.data.gen = 3, n.reps = 10, save.results = FALSE, k.pen = k.pen,
                    seed.metric = TRUE, check.convergence = FALSE, roe = FALSE, fit.metric = "bic", correction = "knp")
}
\arguments{
  \item{input.data}{must be data frame, first column is
  treated as dependent variable.}

  \item{output}{string argument for use in naming file
  output. \code{\link{gemmEst}} writes a .RData file in the
  current working directory each time the function is
  called.}

  \item{n.beta}{Number of beta vectors to generate per
  replication.}

  \item{p.est}{Percept of data used to estimate the model.
  Values less than 1 will cause gemmModel to produce
  cross-validation estimates.}

  \item{n.data.gen}{Number of times the fitting process
  will be repeated.}

  \item{n.reps}{Number of generations for
  \code{\link{genAlg}}.}

  \item{save.results}{Logical value to determine whether
  the resulting gemm object is saved to a .RData file.}

  \item{k.pen}{Penalty term for BIC, as calculated by
  \code{\link{gemm.formula}}.}

  \item{seed.metric}{Logical value to control whether
  \code{\link{genAlg}} is seeded with OLS regression
  weights or random values.}

  \item{check.convergence}{Logical value to indicate
  whether BIC for each generation is retained, mostly
  useful for checking performance of \code{\link{genAlg}}.}

  \item{roe}{Logical value to determine whether region of
  equivalence data are retained.}
  
  \item{fit.metric}{Whether or not to reorder equivalent BIC-tau values based on
  minimized least-squares error.}
  
  \item{correction}{PLACEHOLDER TEXT.}
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
\keyword{ordinal}

