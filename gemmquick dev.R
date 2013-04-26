##### gemmR #####
#
# Authors: Jeff Chrabaszcz & Joe Tidwell
#
# An R implementation of the General Monotone Model, (Dougherty & Thomas, 2012),
# with some improvements.
#
# Please cite:
#
# @article{dougherty2012robust,
#   title={Robust decision making in a nonlinear world},
#   author={Dougherty, Michael R and Thomas, Rick P},
#   journal={Psychological Review},
#   volume={119},
#   number={2},
#   pages={321--344},
#   year={2012},
#   publisher={American Psychological Association},
# }
#
##### TODO #####
# * better plot function (ask argument)
# * error messages
# * GA optimization
# * gemmFit optimization
# * summary
# * predict?
# * categorical predictors (probably involves using factor properly) **Joe**
##### Bugs #####
#
##### Ideas #####
# * posterior predictive checks
################################################################################

##### Dependencies #####

##### GeMM Functions #####

geneticAlgorithm <- function(metric.beta, n.beta, n.super.elites, p, reps,
  bestmodels, seed.metric) {
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
#   seed.metric    - control whether lm() estimated betas seed GA. Default is  #
#                    TRUE                                                      #
################################################################################
# ls would control whether to seed the betas with OLS estimates. Currently this
# happens by default and cannot be changed outside of the function.
  if (seed.metric != TRUE) {metric.beta <- runif(p)}
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
# turn half of all candidate betas negative
#   temp.rand <- rbinom(prod(dim(betas)), 1, .5)
#   temp.rand[temp.rand == 0] <- -1
#   temp.rand <- matrix(temp.rand, ncol = ncol(betas))
#   betas <- betas * temp.rand
  y <- betas[1:n.beta,]
  return(y)
}

gemmFit <- function(n, betas, data, p, k.cor, pearson) {
################################################################################
# Function generates model estimates based on sets of weights and predictors,  #
# calculates Kendall's tau between dependent variable and model predictions.   #
#   n     - number of betas for which fit statistics will be calculated.       #
#   betas - matrix of betas, rows are different collections of betas, columns  #
#           are different predictors.                                          #
#   data  - original predictors and outcome used to calculate fit.             #
#   p     - number of predictors.                                              #
#                                                                              #
#  NOTE: k or p for tau transformation needs to be solved.                     #
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
    if (pearson) {
      r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
    }
  }
# this might cause problems, reverses the scale for any negative correlations
# and recalculates fit. Might be able to just multiply by -1?
  if (tau < 0) {
    betas <- betas * -1
    tau <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)),
      method = "kendall")
    if (pearson) {
      r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
    }
  }
  knp <- sin(pi/2*tau*((n-p-1)/n))
  bic <- n * log(1 - knp ^ 2) + k.cor * log(n)
  y <- list(bic = bic)
  if (pearson) {
    y$r <- r
    y$tau <- tau
  }
  return(y)
}

gemmEst <- function(input.data, output = "gemmr", n.beta = 2000, p.est = 1,
  n.data.gen = 3, n.reps = 10, save.results = FALSE, k.pen = k.pen,
  seed.metric = TRUE, check.convergence = FALSE) {
################################################################################
# Function controls the GeMM process. Takes data and, over successive          #
# replications, uses geneticAlgorithm to generate candidate beta vectors,      #
# calculates ordinal model fit using these betas, and produces an output that  #
# reports weights and fit statistics for best models at each generation,       #
# (optionally, for cross-validation as well).                                  #
#   input.data  - must be data frame, first column is treated as dependent     #
#                 variable.                                                    #
#   output      - string argument for use in naming file output. gemmModel     #
#                 writes a .RData file in the current working directory each   #
#                 time the function is called.                                 #
#   n.beta      - Number of beta vectors to generate per replication. Default  #
#                 is 2000.                                                     #
#   p.est       - Percept of data used to estimate the model. Default is 1,    #
#                 values less than 1 will cause gemmModel to produce           #
#                 cross-validation estimates.                                  #
#   n.data.gen  - Number of times the entire GeMM process will be repeated,    #
#                 due for removal.                                             #
#   n.reps      - Number of replications, default is 10.                       #
#   k.pen       - additional penalty to BIC for including, NA by default.      # 
#   seed.metric - control whether lm() estimated betas seed GA. Default is     #
#                 TRUE                                                         #
################################################################################
  bestmodels <- c()
  var.name <- colnames(input.data[,-1])
  n.super.elites <- round(n.beta/16)
  fit.out <- matrix(rep(0, times = n.data.gen * (dim(input.data)[2])),
    nrow = n.data.gen)
  fit.out.r <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
  fit.out.tau <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
  if (check.convergence) {
    converge.bic <- matrix(rep(0, times = (n.reps * n.data.gen)),
      ncol = n.data.gen)
    converge.beta <- matrix(rep(0,
      times = (n.reps * n.data.gen * (dim(input.data)[2] - 1))),
      ncol = (dim(input.data)[2] - 1))
  }
  if (p.est < 1) {
    gemm.cross.out <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
    gemm.cross.out.r <- gemm.cross.out
    gemm.cross.out.tau <- gemm.cross.out
  }
  for (datagen in 1:n.data.gen) {
    get.r <- FALSE
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
      if (reps == n.reps) {
        get.r <- TRUE
      }
# beta generation here
      betas <- geneticAlgorithm(metric.beta, n.beta, n.super.elites, p, reps,
        bestmodels, seed.metric)
      betas <- as.matrix(betas)
# calculate penalized k for interactions
      k.cor <- rep(1, times = nrow(betas))
      if (!is.null(dim(k.pen))) {
        k.cor <- kCorFact(k.pen, betas)
        k.cor <- matrix(k.cor, ncol = 1)
      }
      fit.stats <- matrix(rep(0, times = (dim(betas)[1])), ncol = 1)
      if (get.r) {
        fit.stats.r <- fit.stats
        fit.stats.tau <- fit.stats
      }
# this loop calculates fit. Could be optimized.
      for (i in 1:dim(betas)[1]) {
        gemm.fit.out <- gemmFit(n, betas[i,], data, p, k.cor[i], pearson = get.r)
        fit.stats[i,] <- gemm.fit.out$bic
        if (get.r) {
          fit.stats.r[i,] <- gemm.fit.out$r
          fit.stats.tau[i,] <- gemm.fit.out$tau
        }
      }
      model.stats <- cbind(fit.stats, betas)
      model.stats <- rbind(rep(0, times = length(model.stats[1,])), model.stats)
      model.stats <- model.stats[order(model.stats[,1]),]
      if (get.r) {
        fit.stats.r <- fit.stats.r[order(model.stats[,1])]
        fit.stats.tau <- fit.stats.tau[order(model.stats[,1])]
      }
      bestmodels <- model.stats[1:(4*n.super.elites),]
      if (check.convergence) {
        converge.bic[reps, datagen] <- bestmodels[1,1]
        converge.beta[(reps + reps * (datagen - 1)),] <- bestmodels[1,-1]
      }
    }
    fit.out[datagen,] <- bestmodels[1,]
    fit.out.r[datagen,] <- fit.stats.r[1]
    fit.out.tau[datagen,] <- fit.stats.tau[1]
    if (p.est < 1) {
      temp.out <- gemmFit(n, betas[i,], cross.val, p, pearson = get.r)
      gemm.cross.out[datagen,] <- temp.out$bic
      gemm.cross.out.r[datagen,] <- temp.out$r
      gemm.cross.out.tau[datagen,] <- temp.out$tau
    }
  }
  coefficients <- matrix(fit.out[,-1], ncol = p,
    dimnames = list(c(), c(colnames(input.data))[-1]))
  fitted.values <- coefficients[as.numeric((fit.out[,1] == 
      min(fit.out[,1]))[1])] * input.data[,-1]
  sim.results <- list(date = date(),
                      call = match.call(),
                      coefficients = coefficients,
                      fitted.values = fitted.values,
                      residuals = unlist(input.data[1] - fitted.values),
                      rank.residuals = rank(input.data[1] -
                          rank(fitted.values)),
                      est.bic = fit.out[,1],
                      est.r = c(fit.out.r),
                      est.tau = c(fit.out.tau),
                      metric.betas = metric.beta,
                      p.vals = p.vals,
                      model = data.frame(input.data))
  if (p.est < 1) {
    sim.results$cross.val.bic <- c(gemm.cross.out)
    sim.results$cross.val.r <- c(gemm.cross.out.r)
    sim.results.cross.val.tau <- c(gemm.cross.out.tau)
  }
  if (check.convergence) {
    sim.results$converge.bic <- converge.bic
    sim.results$converge.beta <- converge.beta
    attr(sim.results, "converge.check") <- TRUE
  }
  if (save.results) {
    save(sim.results, file = paste(output, ".Rdata"))
  }
  return(sim.results)
}

kCorFact <- function(k.pen, beta.vecs) {
################################################################################
# Computes the correctd k for interaction terms. Returns a vector.             #
#   k.pen     - factor matrix from model.frame called in gemm.formula.         #
#   beta.vecs - matrix of beta vectors produced by geneticAlgorithm.           #
################################################################################
  k <- dim(k.pen)[2]
  factors <- log2(k+1)
  levels <- rbind(k.pen, diag(k)[(factors+1):k,])
  return(apply(beta.vecs,1, function(x) sum(levels%*%x!=0)))
}
################################################################################

##### Package functions #####
gemm <- function(x, ...) UseMethod("gemm")

gemm.default <- function(x, k.pen = k.pen, ...) {
  est <- gemmEst(input.data = x, k.pen = k.pen, ...)
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
  est <- gemm.default(cbind(y, x), k.pen = k.pen, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

plot.gemm <- function(x, ...) {
  par(mfrow = c(1,2))
  if (!is.null(attr(x, "converge.check"))) {
    par(mfrow = c(1,3))
    convergencePlot(x$converge.bic)
  }
  plot(rank(fitted.values(x)), rank(x$model[1]))
  plot(fitted.values(x), unlist(x$model[1]))
}

convergencePlot <- function(beta, ...) {
  chains <- ncol(beta)
  max.rep <- nrow(beta)
  xrange <- c(1, max.rep) 
  yrange <- c(min(beta), max(beta))
  plot(xrange, yrange, type="n", xlab = "rep #", ylab= "BIC")
  colors <- rainbow(chains) 
  linetype <- c(1:chains) 
  plotchar <- seq(1:chains)
  for (i in 1:chains) {
    lines(1:max.rep, beta[,i], type="b", lwd=1.5,
      lty=linetype[i], col=colors[i], pch=plotchar[i]) 
  }  
  title("Convergence of BICs")
  legend(xrange[1], yrange[2], 1:chains, cex=0.8, col=colors, pch=plotchar,
    lty=linetype, title="Chains")
}