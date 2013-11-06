\name{genAlg}
\alias{genAlg}
\title{genAlg}
\usage{
  genAlg(metricbeta, nbeta, nsuperelites, p, reps,
    bestmodels, seedmetric)
}
\arguments{
  \item{metricbeta}{Weights derived using multiple
  regression. Overwritten if \code{seedmetric} is FALSE.}

  \item{nbeta}{Number of candidate weight vectors in each
  generation.}

  \item{nsuperelites}{Number of candidate weight vectors to
  involve in permutation for reps > 1.}

  \item{p}{Number of potential predictors.}

  \item{reps}{Generation number. For reps == 1, entirely
  new weights are generated. When reps > 1, bestmodels are
  used to generate permutations.}

  \item{bestmodels}{Matrix of best candidate weight vectors
  from previous generation.}

  \item{seedmetric}{If TRUE, multiple regression weights
  are used to seed the genetic algorithm. Otherwise, random
  weights are used.}
}
\value{
  Returns a matrix of candidate weights with rows for each predictor and columns
  for each unique vector of betas.
}
\description{
  Genetic algorithm for gemmR. Currently has fixed scaling factors so predictors
  should be normalized. Heavily seeded with zero values to interact properly
  with AIC/BIC calculation.
}
\examples{
}

