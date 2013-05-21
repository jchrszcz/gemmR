gemmR
=====

General Monotone Model in R. This is an R version of the non-parametric regression algorithm proposed by Dougherty & Thomas (2012). The algorithm works as follows:

1. Generate random weights with a genetic algorithm
2. Test the fit between the sum of the weighted predictors and an outcome using Bayesian Information Criterion based on tau-to-r transformation.
3. Select top 1/16th of candidate beta vectors, mutate and generate new weights, repeat for some number of generations.
4. Repeat entire process some number of times to check for convergence.

The resulting `gemm` object is currently an S3 class object. We're working on reducing runtime, (mostly through parallelization at this point), and writing a proper helpfile. Information on the basis for GeMM can be found in the [original paper](http://www.bsos.umd.edu/psyc/dougherty/pdf%20articles/DoughertyThomas2012Rev.pdf).
