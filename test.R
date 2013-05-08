require(Rcpp)
sourceCpp("new_ga.cpp")
metricbeta <- matrix(c(1,1,1), ncol = 3)
nbeta <- 100
nsuperelites <- 6
p <- 3
reps <- 2
bestmodels <- matrix(rnorm(96), nrow = 4)
genAlg(metricbeta, nbeta = 2000, nsuperelites, p, reps, bestmodels, seedmetric = TRUE)

cppFunction('
  NumericVector randTest() {
    NumericVector parent1(1);
    NumericVector parent2(1);
    int nsuperelites = 6;
    parent1 = floor(nsuperelites * runif(1) + .5);
    do { (parent2 = floor(nsuperelites * runif(1) + .5));
      } while (parent1(0) == parent2(0));
    return parent1;
    return parent2;
  }
')

##### Testing

metric.beta <- t(metricbeta)
best.models <- t(bestmodels)

system.time(
  for (i in 1:1000) {
    genAlg(metricbeta, nbeta = 2000, nsuperelites, p, reps, bestmodels, seedmetric = FALSE)
  }
)

system.time(
  for (i in 1:1000) {
    geneticAlgorithm(metric.beta, n.beta = 2000, n.super.elites=nsuperelites, p, reps, best.models, seed.metric = FALSE)
  }
)
