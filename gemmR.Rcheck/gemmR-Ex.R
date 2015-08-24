pkgname <- "gemmR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('gemmR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("convergenceplot")
### * convergenceplot

flush(stderr()); flush(stdout())

### Name: convergencePlot
### Title: convergencePlot function for visualizing genetic algorithm in
###   gemmR
### Aliases: convergencePlot

### ** Examples

  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + hp, data = mtcars, check.convergence = TRUE,
    seed.metric = FALSE, n.chains = 3, n.gens = 3, n.beta = 200)
  with(gemm.model, convergencePlot(converge.fit.metric, fit.metric))



cleanEx()
nameEx("gemm")
### * gemm

flush(stderr()); flush(stdout())

### Name: gemm
### Title: Fit General Monotone Models.
### Aliases: gemm gemm.default gemm.formula print.gemm plot.gemm
###   summary.gemm print.summary.gemm nobs.gemm deviance.gemm logLik.gemm
### Keywords: ordinal, regression

### ** Examples

  ## Not run: 
##D     data(culture)
##D     gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars, check.convergence = TRUE)
##D     print(gemm.model)
##D     plot(gemm.model)
##D   
## End(Not run)



cleanEx()
nameEx("gemmEst")
### * gemmEst

flush(stderr()); flush(stdout())

### Name: gemmEst
### Title: Fit General Monotone Models.
### Aliases: gemmEst
### Keywords: ordinal, regression

### ** Examples

  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars,
    check.convergence = TRUE, n.chains = 3, n.gens = 3, n.beta = 200)
  print(gemm.model)
  plot(gemm.model)



cleanEx()
nameEx("genAlg")
### * genAlg

flush(stderr()); flush(stdout())

### Name: genAlg
### Title: genAlg
### Aliases: genAlg

### ** Examples

  p <- 4
  gen.alg <- genAlg(matrix(rnorm(p), nrow = p), 5, 2, p, 1, matrix(1), TRUE)



cleanEx()
nameEx("list2gemm")
### * list2gemm

flush(stderr()); flush(stdout())

### Name: list2gemm
### Title: Collapse a uniform list of 'gemm' objects.
### Aliases: list2gemm
### Keywords: ordinal, regression

### ** Examples

  library(gemmR)
  data(mtcars)
  fit <- list()
  for (i in 1:3) {
      fit[[i]] <- gemm(mpg ~ disp + cyl, data = mtcars, n.chains = 1, n.gens = 3, n.beta = 200)
  }
  gemm.model <- list2gemm(fit)
  summary(gemm.model)



cleanEx()
nameEx("predict.gemm")
### * predict.gemm

flush(stderr()); flush(stdout())

### Name: predict.gemm
### Title: Predict method for general monotone models.
### Aliases: predict.gemm
### Keywords: ordinal, regression

### ** Examples

  data(mtcars)
  gemm.model <- gemm(mpg ~ disp + cyl, data = mtcars, check.convergence = TRUE, n.beta = 200)
  predict(gemm.model, tie.struct = TRUE)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
