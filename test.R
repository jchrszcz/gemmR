metricbeta <- matrix(c(1,1,1), ncol = 3)
nbeta <- 10
nsuperelites <- 2
p <- 3
reps <- 1
bestmodels <- matrix(c(-1,1,1,1), ncol = 4)
genAlg(metricbeta, nbeta, nsuperelites, p, reps = 1, bestmodels, seedmetric = FALSE)