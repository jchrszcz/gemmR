metricbeta <- matrix(c(1,1,1), ncol = 3)
nbeta <- 10
nsuperelites <- 2
p <- 3
reps <- 1
bestmodels <- matrix(rep(1, times = 1000), nrow = 4)
genAlg(metricbeta, nbeta = 200, nsuperelites, p, reps = 2, bestmodels, seedmetric = FALSE)