#TODO b <- as.numeric(y!=0) for dichotomization may be much faster
# all crossval runs
# search all possible models in analytic space


dichotomizeData <- function(x) {
  if (is.data.frame(x) != T) {stop("input must be a dataframe")}
  y <- as.data.frame(matrix(rep(0, times = (dim(x)[1]*dim(x)[2])),
                            ncol = length(x)))
  for (i in 1:length(x)) {
    y[,i] <- ifelse(x[,i] > median(x[,i]), 1, 0)
    y[,i] <- y[,i] + rnorm(length(y[,i])) * sqrt(10E-6)
  }
  names(y) <- names(x)
  return(y)
}

normalizeData <- function(x, rank) {
  # possible values for "rank" in normalizeData are:
  # c("orig","norm","rank","nlog","sqrt")
  if (is.data.frame(x) != T) {stop("input must be a dataframe")}
  size <- dim(x)
  if (rank == "norm") {
    for (i in 1:size[2]) {
      x[,i] <- (x[,i] - min(x[,i]))/(max(x[,i]) - min(x[,i]))
    }
  }
  if (rank == "rank") {
    for (i in 1:size[2]) {
      x[,i] <- rank(x[,i], ties.method = "average")
    }
  }
  if (rank == "nlog") {
    if (sum(sum(x < 0)) > 0)
      stop("negative numbers exist in untransformed data")
    for (i in 1:size[2]) {
      x[,i] <- log(x[,i] + 1)
    }
  }
  if (rank == "sqrt") {
    for (i in 1:size[2]) {
      x[,i] <- sqrt(x[,i] + 1)
    }
  }
  return(x)
}

geneticAlgorithm <- function(metric.beta, n.beta, n.beta.vectors,
                             n.super.elites, p, reps, weights, restrict.parameters) {
  ls <- 1
  if (ls != 1) {metric.beta <- rnorm(p)}
  if (reps == 1) {
    betas <- matrix(rep(0, times = n.beta * p), ncol = p)
    scaling <- sqrt(.1)
    betas[1,] <- metric.beta
    for (i in 2:n.beta) {
      temp.rand <- runif(p)
      temp.norm <- rnorm(p)
      if (i >= 2 & i < 1000) {
        betas[i,] <- ifelse(temp.rand < .5, 1, 0)
        betas[i,] <- betas[i,] * metric.beta # check this
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
                                       1 + temp.norm.2[,2] * scaling, -1 + temp.norm.2[,2] * scaling)
      }
    }
  }
  if (reps > 1) {                                # genetic algorithm starts here
    load("bestmodels.Rdata")
    size <- dim(bestmodels)
    elites <- bestmodels
    sorted.elites <- elites[order(elites[,1]),]
    super.elites <- sorted.elites[1:n.super.elites,]
    #     cumulative.fit <- rep(0, times = size[1])
    #     sum.fit <- sum(abs(sorted.elites[,1]))
    #     cumulative.fit[1] <- abs(sorted.elites[1,1])/sum.fit
    #     for (i in 2:size[1]) {
    #       cumulative.fit[i] <- cumulative.fit[i-1] + abs(sorted.elites[i,1])/sum.fit
    #     }
    #    # next bit could be optimized
    #     temp.betas.a <- c()
    #     temp.betas.b <- c()
    #     for (i in 1:size[1]) {
    #       p.rand <- runif(1)
    #       for (j in 2:size[1]) {
    #         if (p.rand >= cumulative.fit[j-1] & p.rand < cumulative.fit[j]) {
    #           temp.betas.a[i,] <- sorted.elites[j,(size[2] - p + 1):size[2]]
    #         }
    #         if (p.rand == 1) {
    #           temp.betas.a[i,] <- sorted.elites[size[1],(size[2] - p + 1):size[2]]
    #         }
    #       }
    #       p.rand <- runif(1)
    #       for (j in 2:size[1]) {
    #         if (p.rand >= cumulative.fit[j-1] & p.rand < cumulative.fit[j]) {
    #           temp.betas.b[i,] <- sorted.elites[j,(size[2] - p + 1):size[2]]
    #         }
    #         if (p.rand == 1) {
    #            temp.betas.b[i,] <- sorted.elites[size[1],(size[2] - p + 1):size[2]]
    #          }
    #        }
    #      }
    temp.betas.a <- sorted.elites[,(size[2] - p + 1):size[2]]
    temp.betas.b <- sorted.elites[,(size[2] - p + 1):size[2]]
    parent.1 <- round(1 + (size[1] - 1) * runif(1))
    parent.2 <- 0
    while (parent.1 == parent.2 | parent.2 == 0) {
      parent.2 <- round(1 + (size[1] - 1) * runif(1))
    }    # Gene crossover (line 171)
    temp.rand <- runif(n.beta/2)
    new.X1 <- matrix(rep(0, times = (n.beta/2 * p)), ncol = p)
    new.X2 <- new.X1
    for (i in 1:(n.beta/2)) { # this loop could be optimized
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
          new.X1[i,] <- as.numeric(c(temp.betas.a[parent.1,1:k], temp.betas.b[parent.2, ((k+1):p)]))
          new.X2[i,] <- as.numeric(c(temp.betas.b[parent.1,1:k], temp.betas.a[parent.2, ((k+1):p)]))
        }
      }
      if (temp.rand[i] >= .85) {
        new.X1[i,] <- as.numeric(temp.betas.a[parent.1,])
        new.X2[i,] <- as.numeric(temp.betas.b[parent.2,])
      }
    }
    temp.rand <- matrix(runif(n.beta*p), ncol = (p*2)) # this should be fixed, bias in X1/X2 selection
    temp.rand.2 <- matrix(runif(n.beta*p), ncol = (p*2))
    for (i in 1:p) {
      new.X1[,i] <- ifelse(temp.rand[,i] < .01, temp.rand[,(i + p)], new.X1[,i])
      new.X2[,i] <- ifelse(temp.rand.2[,3] < .01, temp.rand.2[,(i + p)], new.X2[,i])
    }
    #new.X3 <- new.X1         # Best Offspring ##### Is this functioning in current code?
    #new.X4 <- new.X2         # need to add X3 and X4 once question is answered
    #
    super.elites <- super.elites[,(7:(p+6))]
    betas <- rbind(as.matrix(super.elites), new.X1, new.X2)
    #if (reps > 1) {
    #  betas <- rbind(super.elites, new.X3, new.X4)
    #}
  }
  #n.beta.vectors <- dim(betas)[1]
  #ncol<- dim(betas)[2]
  #temp.betas <- betas
  #for (i in 1:ncol) {
  #  temp.betas[,i] <- ifelse(temp.betas[,i] != 0, 1, 0)
  #}
  # skipping restrictparameters and weights for now
  y <- betas[1:n.beta,]
  return(y)
}

gemmFit <- function(n, betas, predictors, criterion, fit.criterion,
                    model.selection, p) {
  fit.criterion <- tolower(fit.criterion)
  model.selection <- tolower(model.selection)
  u.predictors <- rep(0, times = n)
  for (i in 1:n) {
    u.predictors[i] <- sum(betas * predictors[i,])
  }
  k <- sum(betas != 0)
  if (cor(criterion, u.predictors) < 0 & is.na(cor(criterion, u.predictors)) == F) {
    betas <- betas * -1
    for (i in 1:n) {
      u.predictors[i] <- sum(betas * predictors[i,])
    }
  }
  if (k == 0) {
    tau <- 0
  }
  ken <- sin(pi/2*tau)
  knp <- sin(pi/2*tau*((n-p-1)/n))
  knk <- sin(pi/2*tau*((n-k-1)/n))
  if (sd(u.predictors) != 0) { 
    pea <- cor(criterion, u.predictors)
  }
  if (sd(u.predictors) == 0) { 
    pea <- 0
  }
  #if (is.na(pea)) {pea <- 0}
  pc1 <- con/(con + dis + ties)
  gam <- (con - dis)/(con + dis)
  success <- (con + .5 * ties)/(con + dis + ties)
  kn2 <- 1 - ((1 - ken^2)*(n-1)/(n - k - 1))
  if (fit.criterion == "pc1" & model.selection == "bc1") {
    bic <- n * log(1 - pc1) + k * log(n)
  }
  if (fit.criterion == "pc1" & model.selection == "sat") {
    bic <- 1 - pc1
  }
  if (fit.criterion == "pc1" & model.selection == "aic") {
    bic <- n * log(1 - pc1) + 2 * k
  }
  if (fit.criterion == "ols" & model.selection == "bc1") {
    bic <- n * log(1 - pea ^ 2) + k * log(n)
  }
  if (fit.criterion == "ols" & model.selection == "sat") {
    bic <- 1 - pea ^ 2
  }
  if (fit.criterion == "ols" & model.selection == "aic") {
    bic <- n * log(1 - pea ^ 2) + 2 * log(n)
  }
  if (fit.criterion == "ols" & model.selection == "ric") {
    bic <- n * log(1 - pea ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "ken" & model.selection == "bc1") {
    bic <- n * log(1 - ken ^ 2) + k * log(n)
  }
  if (fit.criterion == "ken" & model.selection == "sat") {
    bic <- 1 - ken
  }
  if (fit.criterion == "ken" & model.selection == "ric") {
    bic <- n * log(1 - ken ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "ken" & model.selection == "aic") {
    bic <- n * log(1 - ken ^ 2) + 2 * k
  }
  if (fit.criterion == "knp" & model.selection == "bc1") {
    bic <- n * log(1 - knp ^ 2) + k * log(n)
  }
  if (fit.criterion == "knp" & model.selection == "sat") {
    bic <- 1 - knp ^ 2
  }
  if (fit.criterion == "knp" & model.selection == "ric") {
    bic <- n * log(1 - knp ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "knp" & model.selection == "aic") {
    bic <- n * log(1 - knp ^ 2) + 2 * k
  }
  if (fit.criterion == "knk" & model.selection == "bc1") {
    bic <- n * log(1 - knk ^ 2) + k * log(n)
  }
  if (fit.criterion == "knk" & model.selection == "sat") {
    bic <- 1 - knk ^ 2
  }
  if (fit.criterion == "knk" & model.selection == "ric") {
    bic <- n * log(1 - knk ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "knk" & model.selection == "aic") {
    bic <- n * log(1 - knk ^ 2) + 2 * k
  }
  if (fit.criterion == "kn2" & model.selection == "bc1") {
    bic <- n * log(1 - kn2) + k * log(n)
  }
  if (fit.criterion == "kn2" & model.selection == "sat") {
    bic <- 1 - kn2
  }
  if (fit.criterion == "kn2" & model.selection == "ric") {
    bic <- n * log(1 - kn2) + 2 * k * log(p)
  }
  if (fit.criterion == "kn2" & model.selection == "aic") {
    bic <- n * log(1 - kn2) + 2 * k
  }
  raw.data <- c(con, dis, ties, crit, both)
  y <- list(fit = c(bic, pc1, success, tau, pea, k), 
            raw.data = c(con, dis, ties, crit, both))
  return(y)
}

gemmCross <- function(n, k, criterion, u.predictors, fit.criterion,
                      model.selection, p) {
  dyads <- (n * (n - 1))/2
  nobs <- n - 1
  con <- 0
  dis <- 0
  both <- 0
  cue <- 0
  crit <- 0
  for (i in 1:nobs) {
    for (j in i:nobs) {
      if (criterion[i] < criterion[j + 1] & u.predictors[i] < u.predictors[j + 1]) {
        con <- con + 1
      }
      if (criterion[i] > criterion[j + 1] & u.predictors[i] > u.predictors[j + 1]) {
        con <- con + 1
      }
      if (criterion[i] < criterion[j + 1] & u.predictors[i] > u.predictors[j + 1]) {
        dis <- dis + 1
      }
      if (criterion[i] > criterion[j + 1] & u.predictors[i] < u.predictors[j + 1]) {
        dis <- dis + 1
      }
      if (criterion[i] != criterion[j + 1] & u.predictors[i] == u.predictors[j + 1]) {
        cue <- cue + 1
      }
      if (criterion[i] == criterion[j + 1] & u.predictors[i] != u.predictors[j + 1]) {
        crit <- crit + 1
      }
      if (criterion[i] == criterion[j + 1] & u.predictors[i] == u.predictors[j + 1]) {
        both <- both + 1
      }
    }
  }
  ties <- cue + crit + both
  tau <- (con - dis)/dyads
  if (k == 0) {
    tau <- 0
  }
  ken <- sin(pi/2*tau)
  knp <- sin(pi/2*tau*((n-p-1)/n))
  knk <- sin(pi/2*tau*((n-k-1)/n))
  if (sd(u.predictors) != 0) { 
    pea <- cor(criterion, u.predictors)
  }
  if (sd(u.predictors) == 0) { 
    pea <- 0
  }
  pc1 <- con/(con + dis + ties)
  gam <- (con - dis)/(con + dis)
  success <- (con + .5 * ties)/(con + dis + ties)
  kn2 <- 1 - ((1 - ken^2)*(n-1)/(n - k - 1))
  if (fit.criterion == "pc1" & model.selection == "bc1") {
    bic <- n * log(1 - pc1) + k * log(n)
  }
  if (fit.criterion == "pc1" & model.selection == "sat") {
    bic <- 1 - pc1
  }
  if (fit.criterion == "pc1" & model.selection == "aic") {
    bic <- n * log(1 - pc1) + 2 * k
  }
  if (fit.criterion == "ols" & model.selection == "bc1") {
    bic <- n * log(1 - pea ^ 2) + k * log(n)
  }
  if (fit.criterion == "ols" & model.selection == "sat") {
    bic <- 1 - pea ^ 2
  }
  if (fit.criterion == "ols" & model.selection == "aic") {
    bic <- n * log(1 - pea ^ 2) + 2 * log(n)
  }
  if (fit.criterion == "ols" & model.selection == "ric") {
    bic <- n * log(1 - pea ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "ken" & model.selection == "bc1") {
    bic <- n * log(1 - ken ^ 2) + k * log(n)
  }
  if (fit.criterion == "ken" & model.selection == "sat") {
    bic <- 1 - ken
  }
  if (fit.criterion == "ken" & model.selection == "ric") {
    bic <- n * log(1 - ken ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "ken" & model.selection == "aic") {
    bic <- n * log(1 - ken ^ 2) + 2 * k
  }
  if (fit.criterion == "knp" & model.selection == "bc1") {
    bic <- n * log(1 - knp ^ 2) + k * log(n)
  }
  if (fit.criterion == "knp" & model.selection == "sat") {
    bic <- 1 - knp ^ 2
  }
  if (fit.criterion == "knp" & model.selection == "ric") {
    bic <- n * log(1 - knp ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "knp" & model.selection == "aic") {
    bic <- n * log(1 - knp ^ 2) + 2 * k
  }
  if (fit.criterion == "knk" & model.selection == "bc1") {
    bic <- n * log(1 - knk ^ 2) + k * log(n)
  }
  if (fit.criterion == "knk" & model.selection == "sat") {
    bic <- 1 - knk ^ 2
  }
  if (fit.criterion == "knk" & model.selection == "ric") {
    bic <- n * log(1 - knk ^ 2) + 2 * k * log(p)
  }
  if (fit.criterion == "knk" & model.selection == "aic") {
    bic <- n * log(1 - knk ^ 2) + 2 * k
  }
  if (fit.criterion == "kn2" & model.selection == "bc1") {
    bic <- n * log(1 - kn2) + k * log(n)
  }
  if (fit.criterion == "kn2" & model.selection == "sat") {
    bic <- 1 - kn2
  }
  if (fit.criterion == "kn2" & model.selection == "ric") {
    bic <- n * log(1 - kn2) + 2 * k * log(p)
  }
  if (fit.criterion == "kn2" & model.selection == "aic") {
    bic <- n * log(1 - kn2) + 2 * k
  }
  y <- c(bic, pc1, success, tau, pea, k)
  return(y)
}

giveDeciles <- function(x, y) {
  if (length(x) < 10) {
    stop("fewer than 10 values for giveDeciles")
  }
  o <- order(x)
  data <- cbind(x[o], y[o])
  size <- dim(data)
  quant.x <- quantile(data[,1], probs = seq(0,1,.1), type = 8)
  quant.y <- quantile(data[,2], probs = seq(0,1,.1), type = 8)
  dec.mat <- matrix(rep(0, times = (size[1] ^ 2)), ncol = size[1])
  q.x <- rep(0, times = length(x))
  q.y <- q.x
  for (i in 1:size[1]) {
    if (data[i,1] <= quant.x[1]) {q.x[i] <- 1}
    if (data[i,1] > quant.x[1] & data[i,1] <= quant.x[2]) {q.x[i] <- 2}
    if (data[i,1] > quant.x[2] & data[i,1] <= quant.x[3]) {q.x[i] <- 3}
    if (data[i,1] > quant.x[3] & data[i,1] <= quant.x[4]) {q.x[i] <- 4}
    if (data[i,1] > quant.x[4] & data[i,1] <= quant.x[5]) {q.x[i] <- 5}
    if (data[i,1] > quant.x[5] & data[i,1] <= quant.x[6]) {q.x[i] <- 6}
    if (data[i,1] > quant.x[6] & data[i,1] <= quant.x[7]) {q.x[i] <- 7}
    if (data[i,1] > quant.x[7] & data[i,1] <= quant.x[8]) {q.x[i] <- 8}
    if (data[i,1] > quant.x[8] & data[i,1] <= quant.x[9]) {q.x[i] <- 9}
    if (data[i,1] > quant.x[9]) {q.x[i] <- 10}
    if (data[i,2] <= quant.y[1]) {q.y[i] <- 1}
    if (data[i,2] > quant.y[1] & data[i,2] <= quant.y[2]) {q.y[i] <- 2}
    if (data[i,2] > quant.y[2] & data[i,2] <= quant.y[3]) {q.y[i] <- 3}
    if (data[i,2] > quant.y[3] & data[i,2] <= quant.y[4]) {q.y[i] <- 4}
    if (data[i,2] > quant.y[4] & data[i,2] <= quant.y[5]) {q.y[i] <- 5}
    if (data[i,2] > quant.y[5] & data[i,2] <= quant.y[6]) {q.y[i] <- 6}
    if (data[i,2] > quant.y[6] & data[i,2] <= quant.y[7]) {q.y[i] <- 7}
    if (data[i,2] > quant.y[7] & data[i,2] <= quant.y[8]) {q.y[i] <- 8}
    if (data[i,2] > quant.y[8] & data[i,2] <= quant.y[9]) {q.y[i] <- 9}
    if (data[i,2] > quant.y[9]) {q.y[i] <- 10}
  }
  data <- cbind(data, q.x, q.y)
  for (i in 1:size[1]) {
    dec.mat[data[i,3], data[i,4]] <- dec.mat[data[i,3], data[i,4]] + 1
  }
  j <- 1
  dec.ind <- rep(0, times = 9)
  for (i in 1:(size[1] - 1)) {
    if (data[i,3] != data[(i + 1),3]) {
      dec.ind[j] = i
      j = j + 1
    }
  }
  y = list(dec.mat = t(dec.mat),
           decileken = rep(0, times = 10))
  return(y)
}

gemmModel <- function(input.data, output, n.beta, model.selection,
                      fit.criterion, p.est, restrict.parameters, scale,
                      dichotomize, n.data.gen, n.reps) {
  # currently working for n.data.gen < 9
  n.beta.vectors <- n.beta/4
  n.super.elites <- n.beta.vectors/4
  fit.out <- matrix(rep(0, times = (n.data.gen * (((dim(data)[2] - 1) * 2) + 6))), nrow = n.data.gen)
  gemm.cross.out <- matrix(rep(0, times = (n.data.gen * 6)), nrow = n.data.gen)
  for (datagen in 1:n.data.gen) {
    data <- input.data
    data <- normalizeData(data, scale)
    if (dichotomize == 1) {
      data <- dichotomizeData(data)
    }
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
    metric.beta <- lin.mod$coef[2:(p + 1)]
    if (datagen == 1) {
      p.vals <- summary(lin.mod)[[4]][-1,4]
    }
    #if (datagen > 1) {
    #  p.vals <- rbind(p.vals, summary(lin.mod)[[4]][-1,4])
    #}
    ps <- ifelse(summary(lin.mod)[[4]][-1,4] < .05, 1, 0)
    for (reps in 1:n.reps) {
      betas <- geneticAlgorithm(metric.beta, n.beta, n.beta.vectors,
                                n.super.elites, p, datagen, 0, 0)
      fit.stats <- matrix(rep(0, times = (dim(betas)[1]*6)), ncol = 6)
      for (i in 1:dim(betas)[1]) {
        gemm.fit.out <- gemmFit(n, betas[i,], data[,-1], data[,1],
                                fit.criterion, model.selection, p)
        fit.stats[i,] <- gemm.fit.out$fit
      }
      model.stats <- cbind(fit.stats, betas)
      model.stats <- rbind(rep(0, times = length(model.stats[1,])), model.stats)
      model.stats <- model.stats[order(model.stats[,1]),]
      bestmodels <- model.stats[1:n.beta.vectors,]
      top.betas <- bestmodels[1,(7:dim(bestmodels)[2])]
      save(bestmodels, file = "bestmodels.Rdata")
      # may need more here
    }
    fit.out[datagen,] <- c(bestmodels[1,], metric.beta)
    if (p.est < 1) {
      cross.val.u.predictors <- rep(0, times = dim(cross.val)[1])
      for (i in 1:length(cross.val.u.predictors)) {
        cross.val.u.predictors[i] <- sum(cross.val[i,-1] * top.betas)
      }
    }
    k <- sum(ifelse(top.betas == 0, 0, 1))
    if (p.est < 1) {
      gemm.cross.out[datagen,] <- gemmCross(dim(cross.val)[1], k, cross.val[,1],
                                            cross.val.u.predictors, fit.criterion, model.selection, (p - 1))
    }
  }
  sim.results <- list(parameters = c(n.beta.vectors, reps, p.est),
                      test.num = c(fit.criterion, model.selection, scale),
                      date = date(),
                      estimation = fit.out,
                      p.vals = p.vals)
  if (p.est < 1) {
    sim.results$cross.val <- gemm.cross.out
  }
  save(sim.results, file = output)
  return(sim.results)
}