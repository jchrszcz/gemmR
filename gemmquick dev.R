##### TODO #####
# * negative weights
# * convergence
# * GA optimization
# * gemmFit optimization
# * S4 class definitions?
# ** summary
# ** predict?
# * correct penalty term for interactions
# * categorical predictors (probably involves using factor properly)
##### Ideas #####
# * force some chains to start without seeding LS estimates to check for
#     robustness to initial conditions?
################################################################################

##### Dependencies #####

##### GeMM Functions #####

geneticAlgorithm <- function(metric.beta, n.beta, n.super.elites, p, reps,
  bestmodels) {
################################################################################
# This functions generates candidate beta weights for each predictor in the    #
# model.                                                                       #
#   metric.beta    - starting beta weights, generally from lm()                #
#   n.beta         - the number of candidate betas for each n.rep. This is     #
#                    controlled by gemmModel()                                 #
#   n.super.elites - argument used to index a certain portion of the beta for  #
#                    different randomization.                                  #
#   p              - number of predictors                                      #
#   reps           - number of replications. If >1, best beta vectors from     #
#                    previous replication are included.                        #
#   bestmodels     - best betas from previous replication.                     #
################################################################################
# ls would control whether to seed the betas with OLS estimates. Currently this
# happens by default and cannot be changed outside of the function.
  ls <- 1
  if (ls != 1) {metric.beta <- rnorm(p)}
# this controls the first generation of beta generation.
  if (reps == 1) {
    betas <- matrix(rep(0, times = n.beta * p), ncol = p)
    scaling <- sqrt(.1)
    betas[1,] <- metric.beta
    for (i in 2:n.beta) {
      temp.rand <- runif(p)
      temp.norm <- rnorm(p)
      if (i >= 2 & i < 1000) {
        betas[i,] <- ifelse(temp.rand < .5, 1, 0)
        betas[i,] <- betas[i,] * metric.beta
      }
      if (i >= 1000 & i < 3000) {
        betas[i,] <- ifelse(temp.rand < .5, 1, betas)
        betas[i,] <- betas[i,] * metric.beta + temp.norm * sqrt(.1)
      }
      if (i >= 3000 & i < 6000) {
        betas[i,] <- ifelse(temp.rand < .5, 1, betas)
        betas[i,] <- betas[i,] * metric.beta + temp.norm * sqrt(.01)
      }
      if (i >= 6000) {
        betas[i,] <- ifelse(temp.rand < .25, 1, ifelse(temp.rand > .75, -1,
                                                       temp.rand))
        betas[i,] <- betas[i,] * metric.beta + temp.norm * sqrt(.5)
      }
    }
    for (i in 1:n.beta) {
      if (sum(betas[i,]) == 0) {
        betas[i,] <- ifelse(betas[i,] < .5, 1, 0)
        temp.norm.2 <- matrix(rnorm(2*length(betas[i,][betas[i,]])), ncol = 2)
        betas[i,][betas[i,]] <- ifelse(temp.norm.2[,1] > 1,
          1 + temp.norm.2[,2] * scaling, -1 +temp.norm.2[,2] * scaling)
      }
    }
  }
# Reps > 1 do not use OLS estimates to seed the model, but do use the top 25%
# (jcz - I think) from previous generation.
  if (reps > 1) {
    size <- dim(bestmodels)
    elites <- bestmodels
    sorted.elites <- elites[order(elites[,1]),]
    super.elites <- sorted.elites[1:n.super.elites,]
    temp.betas.a <- as.matrix(sorted.elites[,2:size[2]])
    temp.betas.b <- as.matrix(sorted.elites[,2:size[2]])
    parent.1 <- round(1 + (size[1] - 1) * runif(1))
    parent.2 <- 0
    while (parent.1 == parent.2 | parent.2 == 0) {
      parent.2 <- round(1 + (size[1] - 1) * runif(1))
    }
    temp.rand <- runif(n.beta/2)
    new.X1 <- matrix(rep(0, times = (n.beta/2 * p)), ncol = p)
    new.X2 <- new.X1
    for (i in 1:(n.beta/2)) {
      parent.1 <- round(1 + (size[1] - 1) * runif(1))
      parent.2 <- 0
      while (parent.1 == parent.2 | parent.2 == 0) {
        parent.2 <- round(1 + (size[1] - 1) * runif(1))
      }
      if (temp.rand[i] < .85) {
        k <- round(1 + ((size[2] - 1) * runif(1)))
        if (k == p) {
          new.X1[i,] <- as.numeric(c(temp.betas.a[parent.1,]))
          new.X2[i,] <- as.numeric(c(temp.betas.b[parent.1,]))
        }
        if (k < p) {
          new.X1[i,] <- as.numeric(c(temp.betas.a[parent.1,1:k],
            temp.betas.b[parent.2, ((k+1):p)]))
          new.X2[i,] <- as.numeric(c(temp.betas.b[parent.1,1:k],
            temp.betas.a[parent.2, ((k+1):p)]))
        }
      }
      if (temp.rand[i] >= .85) {
        new.X1[i,] <- as.numeric(temp.betas.a[parent.1,])
        new.X2[i,] <- as.numeric(temp.betas.b[parent.2,])
      }
    }
    temp.rand <- matrix(runif(n.beta*p), ncol = (p*2))
    temp.rand.2 <- matrix(runif(n.beta*p), ncol = (p*2))
    for (i in 1:p) {
      new.X1[,i] <- ifelse(temp.rand[,i] < .01, temp.rand[,(i + p)], new.X1[,i])
      new.X2[,i] <- ifelse(temp.rand.2[,i] < .01,
        temp.rand.2[,(i + p)], new.X2[,i])
    }
    super.elites <- super.elites[,-1]
    betas <- rbind(as.matrix(super.elites), new.X1, new.X2)
  }
  y <- betas[1:n.beta,]
  return(y)
}

gemmFit <- function(n, betas, data, p) {
################################################################################
# Function generates model estimates based on sets of weights and predictors,  #
# calculates Kendall's tau between dependent variable and model predictions.   #
#   n     - number of betas for which fit statistics will be calculated.       #
#   betas - matrix of betas, rows are different collections of betas, columns  #
#           are different predictors.                                          #
#   data  - original predictors and outcome used to calculate fit.             #
#   p     - number of predictors.                                              #
################################################################################
  
# null models are common, these lines return 0 for both fit metrics when all
# betas are 0.
  if (sum(betas == 0) == p) {
    tau <- 0
    r <- 0
  }
# non-null models trigger these lines of code, which produce the summed products
# of weights and predictors, return the correlation coefficient between those
# predictions and the predicted variable.
  if (sum(betas == 0) != p) {
    tau <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)),
      method = "kendall")
    r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
  }
# this might cause problems, reverses the scale for any negative correlations
# and recalculates fit. Might be able to just multiply by -1?
  if (r < 0) {
    betas <- betas * -1
    tau <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)),
      method = "kendall")
    r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
  }
  k <- sum(betas != 0)
  knp <- sin(pi/2*tau*((n-p-1)/n))
  bic <- n * log(1 - knp ^ 2) + k * log(n)
  y <- list(bic = bic, r = r)
  return(y)
}

gemmEst <- function(input.data, output = "gemmr", n.beta = 2000, p.est = 1,
  n.data.gen = 1, n.reps = 10, save.results = FALSE, k.pen = NA) {
################################################################################
# Function controls the GeMM process. Takes data and, over successive          #
# replications, uses geneticAlgorithm to generate candidate beta vectors,      #
# calculates ordinal model fit using these betas, and produces an output that  #
# reports weights and fit statistics for best models at each generation,       #
# (optionally, for cross-validation as well).                                  #
#   input.data - must be data frame, first column is treated as dependent      #
#                variable.                                                     #
#   output     - string argument for use in naming file output. gemmModel      #
#                writes a .RData file in the current working directory each    #
#                time the function is called.                                  #
#   n.beta     - Number of beta vectors to generate per replication. Default   #
#                is 2000.                                                      #
#   p.est      - Percept of data used to estimate the model. Default is 1,     #
#                values less than 1 will cause gemmModel to produce            #
#                cross-validation estimates.                                   #
#   n.data.gen - Number of times the entire GeMM process will be repeated,     #
#                due for removal.                                              #
#   n.reps     - Number of replications, default is 10.                        #
#   k.pen      - additional penalty to BIC for including, NA by default.       #
################################################################################
  bestmodels <- c()
  var.name <- colnames(input.data[,-1])
  n.super.elites <- round(n.beta/16)
  fit.out <- matrix(rep(0, times = n.data.gen * (dim(input.data)[2])),
    nrow = n.data.gen)
  fit.out.r <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
  if (p.est < 1) {
    gemm.cross.out <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
    gemm.cross.out.r <- gemm.cross.out
  }
  for (datagen in 1:n.data.gen) {
    data <- input.data
    size <- dim(data)
    est.ss <- floor(p.est * size[1])
    data <- as.matrix(data[order(runif(size[1])),])
    if (p.est < 1) {
      cross.val <- data[((est.ss + 1):size[1]),]
    }
    data <- data[1:est.ss,]
    size <- dim(data)
    n <- size[1]
    p <- size[2] - 1
    lin.mod <- lm(data[,1] ~ data[,2:(p + 1)])
    if(sum(is.na(lin.mod))) {
      error("lm() generates NA")
    }
    metric.beta <- lin.mod$coef[2:(p + 1)]
    names(metric.beta) <- names(data[2:length(data)])
    p.vals <- summary(lin.mod)[[4]][-1,4]
    names(p.vals) <- names(data[2:length(data)])
    ps <- ifelse(summary(lin.mod)[[4]][-1,4] < .05, 1, 0)
    for (reps in 1:n.reps) {
# beta generation here
      betas <- geneticAlgorithm(metric.beta, n.beta, n.super.elites, p, reps,
        bestmodels)
      betas <- as.matrix(betas)
      fit.stats <- matrix(rep(0, times = (dim(betas)[1])), ncol = 1)
      fit.stats.r <- fit.stats
# this loop calculates fit. Could be optimized.
      for (i in 1:dim(betas)[1]) {
        gemm.fit.out <- gemmFit(n, betas[i,], data, p)
        fit.stats[i,] <- gemm.fit.out$bic
        fit.stats.r[i,] <- gemm.fit.out$r
      }
      model.stats <- cbind(fit.stats, betas)
      model.stats <- rbind(rep(0, times = length(model.stats[1,])), model.stats)
      model.stats <- model.stats[order(model.stats[,1]),]
      fit.stats.r <- fit.stats.r[order(model.stats[,1])]
      bestmodels <- model.stats[1:(4*n.super.elites),]
      top.betas <- bestmodels[1,-1]
    }
    fit.out[datagen,] <- bestmodels[1,]
    fit.out.r[datagen,] <- fit.stats.r[1]
    k <- sum(top.betas != 0)
    if (p.est < 1) {
      temp.out <- gemmFit(n, betas[i,], cross.val, p)
      gemm.cross.out[datagen,] <- temp.out$bic
      gemm.cross.out.r[datagen,] <- temp.out$r
    }
  }
  sim.results <- list(date = date(),
                      call = match.call(),
                      coefficients = matrix(fit.out[,-1], ncol = p),
                      est.bic = fit.out[,1],
                      est.r = c(fit.out.r),
                      metric.betas = metric.beta,
                      p.vals = p.vals,
                      var.name = var.name)
  colnames(sim.results$coefficients) <- sim.results$var.name
  if (p.est < 1) {
    sim.results$cross.val.bic <- c(gemm.cross.out)
    sim.results$cross.val.r <- c(gemm.cross.out.r)
  }
  if (save.results) {
    save(sim.results, file = paste(output, ".Rdata"))
  }
  return(sim.results)
}
################################################################################

##### Package functions #####
gemm <- function(x, ...) UseMethod("gemm")

gemm.default <- function(x, k.pen = NA, ...) {
  est <- gemmEst(input.data = x, k.pen, ...)
  class(est) <- "gemm"
  est
}

print.gemm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nBIC:\n")
  print(x$est.bic)
}

gemm.formula <- function(formula, data=list(), ...) {
  mf <- model.frame(formula=formula, data=data)
  # removes the intercept column (intercept isn't meaningful)
  if (attributes(attributes(mf)$terms)$intercept == 1) {
    attributes(attributes(mf)$terms)$intercept <- 0
  }
  # retains factor matrix if any interactions are in the model, (necessary for
  #   correctly penalizing BIC)
  if (sum(attributes(attributes(mf)$terms)$order) >= 1) {
    k.pen <- attributes(attributes(mf)$terms)$factor[-1,]
  }
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- gemm.default(cbind(y, x), k.pen, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}