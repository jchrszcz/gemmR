# License:
# This file is part of gemmR. gemmR is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gemmR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gemmR.  If not, see <http://www.gnu.org/licenses/>.

#' gemmR controller function.
#'
#' This function is called by \code{\link{gemm.default}} and runs a handful of subordinate functions to produce a gemm object. Takes data and, over successive replications, uses geneticAlgorithm to generate candidate beta vectors, calculates ordinal model fit using these betas, and produces an output that reports weights and fit statistics for best models at each generation, (optionally, for cross-validation as well).
#' @param input.data must be data frame, first column is treated as dependent variable.
#' @param output string argument for use in naming file output. \code{\link{gemmEst}} writes a .RData file in the current working directory each time the function is called.
#' @param n.beta Number of beta vectors to generate per replication.
#' @param p.est Percept of data used to estimate the model. Values less than 1 will cause gemmModel to produce cross-validation estimates.
#' @param n.data.gen Number of times the fitting process will be repeated.
#' @param n.reps Number of generations for \code{\link{genAlg}}.
#' @param save.results Logical value to determine whether the resulting gemm object is saved to a .RData file.
#' @param k.pen Penalty term for BIC, as calculated by \code{\link{gemm.formula}}.
#' @param seed.metric Logical value to control whether \code{\link{genAlg}} is seeded with OLS regression weights or random values.
#' @param check.convergence Logical value to indicate whether BIC for each generation is retained, mostly useful for checking performance of \code{\link{genAlg}}.
#' @param roe Logical value to determine whether region of equivalence data are retained.
#' @keywords ordinal
#' @export

gemmEst <- function(input.data, output = "gemmr", n.beta = 8000, p.est = 1,
                    n.data.gen = 3, n.reps = 10, save.results = FALSE, k.pen = k.pen,
                    seed.metric = TRUE, check.convergence = FALSE, roe = FALSE) {
  bestmodels <- matrix(rep(0, times = (ncol(input.data) - 1)))
  var.name <- colnames(input.data[,-1])
  n.super.elites <- round(n.beta/16)
  fit.out <- matrix(rep(0, times = n.data.gen * (dim(input.data)[2])),
                    nrow = n.data.gen)
  fit.out.r <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
  fit.out.tau <- matrix(rep(0, times = n.data.gen), nrow = n.data.gen)
  if (roe) {
    roe.mat <- matrix(0, nrow = (n.beta * n.reps * n.data.gen),
      ncol = ncol(input.data))
  }
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
    metricbeta <- matrix(lin.mod$coef[2:(p + 1)], ncol = p)
    names(metricbeta) <- names(data[2:length(data)])
    p.vals <- summary(lin.mod)[[4]][-1,4]
    names(p.vals) <- names(data[2:length(data)])
    ps <- ifelse(summary(lin.mod)[[4]][-1,4] < .05, 1, 0)
    for (reps in 1:n.reps) {
      if (reps == n.reps) {
        get.r <- TRUE
      }
      # beta generation here
      betas <- genAlg(metricbeta, n.beta, n.super.elites, p, reps,
                 t(bestmodels), seed.metric)
      betas <- t(as.matrix(betas))
      # calculate penalized k for interactions
      k.cor <- rep(1, times = nrow(betas))
      if (!is.null(dim(k.pen))) {
        #k.cor <- kCorFact(k.pen, betas)
        k.cor <- apply(betas,1, function(x) sum(as.matrix(k.pen)%*%x!=0))
        k.cor <- matrix(k.cor, ncol = 1)
      } 
      fit.stats <- matrix(rep(0, times = (dim(betas)[1])), ncol = 1)
      if (get.r) {
        fit.stats.r <- fit.stats
        fit.stats.tau <- fit.stats
      }
      fitStats <- gemmFitRcppI(n, betas, data, p, k.cor, get.r)
      fit.stats <- fitStats$bic
      fit.stats.r <- fitStats$r
      fit.stats.tau <- fitStats$tau      
      model.stats <- cbind(fit.stats, betas)
      model.stats <- rbind(rep(0, times = length(model.stats[1,])), model.stats)
      model.stats <- model.stats[order(model.stats[,1]),]
      if (get.r) {
        fit.stats.r <- fit.stats.r[order(model.stats[,1])]
        fit.stats.tau <- fit.stats.tau[order(model.stats[,1])]
      }
      if (roe) {
      	# check this
        roe.mat[1:n.beta + (n.beta * (reps - 1)) +
          (n.beta * n.reps * (datagen - 1)),] <- model.stats[-1,]
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
      tempOut <- gemmFitRcppI(nrow(cross.val), matrix(bestmodels[1,2:(p+1)], nrow = 1), cross.val, p, k.cor, get.r)
      gemm.cross.out[datagen,] <- tempOut$bic
      gemm.cross.out.r[datagen,] <- tempOut$r
      gemm.cross.out.tau[datagen,] <- tempOut$tau
    }
  }
  coefficients <- matrix(fit.out[,-1], ncol = p,
    dimnames = list(c(), c(colnames(input.data))[-1]))
  best.coef <- matrix(fit.out[fit.out[,1] == min(fit.out[,1]), -1], ncol = p)
  if (nrow(best.coef) > 1) {
    best.coef <- best.coef[1,]
  }
  fitted.values <- matrix(input.data[,-1], ncol = p) %*% matrix(best.coef, ncol = 1)
  if (roe) {
  	roe.df <- data.frame(roe.mat)
  	roe.df$beta <- rep(1:n.beta, times = n.reps * n.data.gen)
  	roe.df$reps <- rep(1:n.reps, each = n.beta, times = n.data.gen)
  	roe.df$data.gen <- rep(1:n.data.gen, each = n.beta * n.reps)
  }
  sim.results <- list(date = date(),
    call = match.call(),
    coefficients = coefficients,
    fitted.values = fitted.values,
    residuals = unlist(input.data[,1] - fitted.values),
    rank.residuals = (rank(input.data[,1]) -
        rank(fitted.values)),
    est.bic = fit.out[,1],
    est.r = c(fit.out.r),
    est.tau = c(fit.out.tau),
    metric.betas = metricbeta,
    p.vals = p.vals,
    model = data.frame(input.data))
  if (p.est < 1) {
    sim.results$cross.val.bic <- c(gemm.cross.out)
    sim.results$cross.val.r <- c(gemm.cross.out.r)
    sim.results.cross.val.tau <- c(gemm.cross.out.tau)
    attr(sim.results, "cross.val") <- TRUE
  }
  if (check.convergence) {
    sim.results$converge.bic <- converge.bic
    sim.results$converge.beta <- converge.beta
    attr(sim.results, "converge.check") <- TRUE
  }
  if (roe) {
    sim.results$roe <- roe.df
    attr(sim.results, "roe") <- TRUE
  }
  if (save.results) {
    save(sim.results, file = paste(output, ".Rdata"))
  }
  return(sim.results)
}

#' gemm 
#'
#' This function is called by \code{\link{gemm.default}} and runs a handful of subordinate functions to produce a gemm object. Takes data and, over successive replications, uses geneticAlgorithm to generate candidate beta vectors, calculates ordinal model fit using these betas, and produces an output that reports weights and fit statistics for best models at each generation, (optionally, for cross-validation as well).
#' @param x Formula input for gemm procedure.
#' @param ... Other arguments to be passed to \code{\link{gemm.formula}} and \code{\link{gemmEst}}.
#' @export

gemm <- function(x, ...) UseMethod("gemm")

#' Default method for gemm
#' 
#' This function is the default method for gemm. It calls \code{\link{gemmEst}}, assigns a "gemm" class to the resultant list object, and returns that list.
#' @param x dataframe to be passed to \code{\link{gemmEst}}. First column is predicted variable.
#' @param k.pen Penalty term, necessary for appropriate penalty related to interaction terms.
#' @param ... Other arguments passed to \code{\link{gemmEst}}.
#' @method gemm default
#' @S3method gemm default
#' @export

gemm.default <- function(x, k.pen = k.pen, ...) {
  est <- gemmEst(input.data = x, k.pen = k.pen, ...)
  class(est) <- "gemm"
  est
}

#' Print method for gemm
#' 
#' This function is the print method for gemm. It returns the most-referenced parts of a gem object: formula call, coefficients, and fit indices.
#' @param x gemm object to be referenced.
#' @param ... Other arguments.
#' @method gemm print
#' @S3method gemm print
#' @export

gemm.print <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nBIC:\n")
  print(x$est.bic)
}

#' Function method for gemm
#' 
#' This is the function method for gemm. It calls \code{\link{gemm.default}} after creating the design matrix and calculating the appropriate penalty term
#' @param formula Linear model equivalent to be fit by gemm.
#' @param data Dataframe referred to by formula.
#' @param ... Other arguments to be passed to \code{\link{gemmEst}}.
#' @method gemm formula
#' @S3method gemm formula
#' @export

gemm.formula <- function(formula, data=list(), ...) {
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)[,-1]
  y <- model.response(mf)
  #main effect variables
  me <- attributes(attributes(mf)$terms)$term.labels[attributes(attributes(mf)$terms)$order==1]
  mm <- model.matrix(attr(mf, "terms"), data=mf)
  fmla <- as.formula(paste("y ~ ", paste(me, collapse= "*")))
  #names 
  names <- data.frame(var=names(mf)[-1])
  names$cnt.betas <- apply(names,1,function(x) ifelse(is.factor(mf[,x]),length(levels(mf[,x]))-1,1))
  vars <- apply(names,1,function(x) if(x[2]==1) {x[1]} else {paste(x[1],1:x[2],sep="")} )
  if(dim(names)[1]==1) {
    names = as.vector(vars)
  } else {
    vars <- lapply(vars,function(x)c("",x))
    mat <- t(unique(expand.grid(vars)))
    names <- apply(mat,2,function(x) paste(x, collapse=':'))
    while (length(grep("::",names))>0) {  
          names <- sub("::",':',names)
    }
    names <- sub(':$', '', names)
    names <- sub('^:', '', names)[-1] 
  }
  #names.betas.all <- attr(terms(fmla),"term.labels")
  names.betas.all <- names
  names.betas.in.model <- attr(terms(mf),"term.labels")
  #skip if only 1 predictor
  if(length(names.betas.all)>1) {
    #full combination of variables/factors
    count <- 0
    for (var in me) {
      if(is.factor(mf[,var])) {
        count <- count + length(levels(mf[,var]))-1
      } else {
        count = count +1
      }
    }
    #identity matrix for all combinations
    lst <- list(NULL)
    for(i in 1:count) {
      lst[[i]] <- c(0,1)  
    }
    lst <- t(expand.grid(lst))[,-1]
    lst <- data.frame(lst)
    lst$assign <- attributes(mm)$assign[-1][1:count]
    #remove columns with within factor comparisons
    tmp <- rep(TRUE, times=dim(lst)[2])
    
    for (i in 1 : max(lst$assign)) {
      tmp2 <- colSums(lst[lst$assign==i,])<2
      tmp <- ifelse(tmp==tmp2 & tmp2==TRUE,TRUE,FALSE)
    }
    k.pen <- lst[,tmp==TRUE]
    colnames(k.pen) <- names.betas.all
    # Select Main Effects
    grep.str <- ""
    for (tmp in me) {
      grep.str <- paste(grep.str,"^",tmp,"([0-9]|)$|",sep="")  
    }
    grep.str <- substr(grep.str, 1, nchar(grep.str)-1)
    keep.main <- grepl(grep.str,  names.betas.all)
    keep <- matrix(keep.main,ncol=length(keep.main))
    # Select interactions 
    interactions <- names.betas.in.model[grepl(":",names.betas.in.model)]
    search.terms <- strsplit(interactions,":")
    if (length(search.terms)>0) {
      keep.tmp <- matrix(ncol=dim(k.pen)[2],nrow=length(search.terms))
      keep.tmp[1,] <- F   
      for (i in 1:length(search.terms)) {
        tmp1 <- !!((apply(t(sapply(search.terms[[i]], grepl, colnames(k.pen),
          ignore.case=TRUE)),2,prod)))
        tmp2 <- unlist(lapply(strsplit(colnames(k.pen),":"),
          function(x) length(x) == length(search.terms[[i]])))
        keep.tmp[i,] <- tmp1*tmp2
      }
      keep <- rbind(keep.main,keep.tmp)
    }
    keep <- ifelse(apply(keep,2,sum)>0,T,F)
    k.pen <- k.pen[,keep]
    k.pen <- k.pen[,order(colSums(k.pen))]
    if(dim(k.pen)[1]!=dim(k.pen)[2]) {
      tmp <- diag(dim(k.pen)[2])[(dim(k.pen)[1]+1):dim(k.pen)[2],]
      # stupid R kludge
      if(!is.null(dim(tmp))) {
        colnames(tmp) <- colnames(k.pen)
      }
      k.pen <- rbind(k.pen, tmp)
    }
  } else {
      k.pen <- matrix(1,ncol=1)
      colnames(k.pen) <- names
  }
  est <- gemm.default(cbind(y, x), k.pen = k.pen, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' plot method for gemm
#' 
#' This function is the plot method for a gemm object.
#' @param x "gemm" class object to be used for plotting. Note that plots differ based on whether \code{check.convergence} was called for the initial gemm fitting.
#' @param ... Other arguments.
#' @method gemm plot
#' @S3method gemm plot
#' @export

gemm.plot <- function(x, ...) {
  par(mfrow = c(1,3))
  if (!is.null(attr(x, "converge.check"))) {
    par(mfrow = c(2,2))
    convergencePlot(x$converge.bic)
  }
  plot(rank(fitted.values(x)), rank(x$model[1]),
    main = "Ordinal model predictions", xlab = "Rank predictions",
    ylab = "Rank criterion")
  plot(fitted.values(x), unlist(x$model[1]), main = "Metric model predictions",
    xlab = "Predictions", ylab = "Criterion")
  plot(order(x$model[1]), x$rank.residuals[order(x$model[1])],
    main = "Rank disparity by criterion rank",
    xlab = "Ordered criterion", ylab = "Rank disparity")
}

#' convergencePlot function for visualizing genetic algorithm in gemmR
#' 
#' This function plots convergence for gemm objects that were called with the "check.convergence" argument. Generally doesn't work well with n.data.gen < 10.
#' @param beta Matrix of BIC values for the best model in each generation.
#' @export

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