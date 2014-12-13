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

gemmEst <- function(input.data, output = "gemmr", n.beta = 8000,
                    n.chains = n.chains, n.gens = 10, save.results = FALSE,
                    k.pen = k.pen, seed.metric = TRUE, 
                    check.convergence = FALSE, roe = FALSE, 
                    fit.metric = fit.metric, correction = "knp", 
                    oclo=TRUE, isTauB = FALSE) {
  # Select fitting function
  getFitMetric <- switch(tolower(fit.metric),
                      bic = function(fitStats) {return(fitStats$bic)},
                      aic = function(fitStats) {return(fitStats$aic)},
                      tau = function(fitStats) {return(1-abs(fitStats$tau))},
                      )

  fit.null <- switch(tolower(fit.metric),
                  bic = 0,
                  tau = 1,
                  aic = 0
                  )
  # Allocate Variables 
  bestmodels <- matrix(rep(0, times = (ncol(input.data) - 1)))
  var.name <- colnames(input.data[,-1])
  n.super.elites <- round(n.beta/16)
  fit.out <- matrix(rep(0, times = n.chains * (dim(input.data)[2])),
                    nrow = n.chains)
  fit.out.r <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.tau <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.tau.a <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.tau.b <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.bic <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.aic <- matrix(rep(0, times = n.chains), nrow = n.chains)
  fit.out.tau.par <- matrix(rep(0, times = n.chains*6), nrow = n.chains)
  colnames(fit.out.tau.par) <- c("pairs","ties.1","ties.2","ties.both","dis","con")

  if (roe) {
    roe.mat <- matrix(0, nrow = (n.beta * n.gens * n.chains),
      ncol = (ncol(input.data) + 1))
  }
  if (check.convergence) {
    converge.fit.metric <- matrix(rep(0, times = (n.gens * n.chains)),
                           ncol = n.chains)
    converge.beta <- matrix(rep(0,
                      times = (n.gens * n.chains * (dim(input.data)[2] - 1))),
                            ncol = (dim(input.data)[2] - 1))
    converge.r <- matrix(rep(0, times = (n.gens * n.chains)),
                            ncol = n.chains)
  }
  for (chains in 1:n.chains) {
    data <- input.data
    size <- dim(data)
    est.ss <- floor(size[1])
    data <- as.matrix(data[order(runif(size[1])),])
    data <- data[1:est.ss,]
    size <- dim(data)
    n <- size[1]
    p <- size[2] - 1
    lin.mod <- lm(data[,1] ~ data[,2:(p + 1)])
    if(sum(is.na(lin.mod$coefficients))) {
      warning("lm() generates NA, seeding with random values")
      seed.metric <- FALSE
    }
    metricbeta <- matrix(lin.mod$coef[2:(p + 1)], ncol = p)
    names(metricbeta) <- names(data[2:length(data)])
    p.vals <- summary(lin.mod)[[4]][-1,4]
    names(p.vals) <- names(data[2:length(data)])
    ps <- ifelse(summary(lin.mod)[[4]][-1,4] < .05, 1, 0)
    for (gens in 1:n.gens) {
      betas <- genAlg(metricbeta, n.beta, n.super.elites, p, gens,
                 t(bestmodels), seed.metric)
      betas <- t(as.matrix(betas))
      k.cor <- rep(1, times = nrow(betas))
      if (!is.null(dim(k.pen))) {
        k.cor <- apply(betas,1, function(x) sum(as.matrix(k.pen)%*%x!=0))
        k.cor <- matrix(k.cor, ncol = 1)
      } 
      fitStats <- gemmFitRcppI(n, betas, data, p, k.cor, correction, isTauB)

      # Depends on which tau is used isTauB_
      if(isTauB) {
        fitStats$tau <- fitStats$tau.b
      } else {
        fitStats$tau <- fitStats$tau.a
      }

      fit.stats <- getFitMetric(fitStats)

      fix.tau <- ifelse(fitStats$tau < 0, -1, 1)
      fitStats$r <- fitStats$r * fix.tau
      fitStats$tau <- fitStats$tau * fix.tau
      fitStats$tau.a <- fitStats$tau.a * fix.tau
      fitStats$tau.b <- fitStats$tau.b * fix.tau
      betas <- betas * fix.tau
      model.stats <- cbind(fit.stats, fitStats$r, betas)
      model.stats <- rbind(c(fit.null,rep(0, times = length(model.stats[1,])-1)), model.stats)
      # Order by BIC, then r if oclo=TRUE
      ifelse(oclo, ord <- order((model.stats[,1]),-model.stats[,2]),
                   ord<-order(model.stats[,1])
            )

      model.stats <- model.stats[ord,]
      fit.stats.r <- c(0, fitStats$r)[ord]
      fit.stats.tau <- c(0, fitStats$tau)[ord]
      fit.stats.tau.a <- c(0, fitStats$tau.a)[ord]
      fit.stats.tau.b <- c(0, fitStats$tau.b)[ord]

      fit.stats.tau.par <- cbind(c(0, fitStats$tau.n.pairs)[ord], c(0, fitStats$tau.n.ties.1)[ord],
                                 c(0, fitStats$tau.n.ties.2)[ord], c(0, fitStats$tau.n.ties.both)[ord],
                                 c(0, fitStats$tau.n.dis)[ord], c(0, fitStats$tau.n.con)[ord])
      fit.stats.bic <- c(0, fitStats$bic)[ord]
      fit.stats.aic <- c(0, fitStats$aic)[ord]
      if (roe) {
      	# check this
        roe.mat[1:n.beta + (n.beta * (gens - 1)) +
          (n.beta * n.gens * (chains - 1)),] <- model.stats[-1,]
      }
      bestmodels <- model.stats[1:(4*n.super.elites),-2]
      if (check.convergence) {
        converge.fit.metric[gens, chains] <- bestmodels[1,1]
        converge.beta[(gens + gens * (chains - 1)),] <- bestmodels[1,-1]
        converge.r[gens, chains] <- model.stats[1,2]
      }
    }
    fit.out[chains,] <- bestmodels[1,]
    fit.out.r[chains,] <- fit.stats.r[1]
    fit.out.tau[chains,] <- fit.stats.tau[1]
    fit.out.tau.a[chains,] <- fit.stats.tau.a[1]
    fit.out.tau.b[chains,] <- fit.stats.tau.b[1]
    fit.out.tau.par[chains,] <- fit.stats.tau.par[1,]
    fit.out.bic[chains,] <- fit.stats.bic[1]
    fit.out.aic[chains,] <- fit.stats.aic[1]
  }
  coefficients <- matrix(fit.out[,-1], ncol = p)
  colnames(coefficients) <- colnames(input.data)[-1]
  coefficients[is.na(coefficients)] <- 0

  # Scale coefficients to metric scale
  y.hats <- matrix(input.data[,-1], ncol = p) %*% t(coefficients) 
  scales <- t(apply(y.hats,2,function(x) coef(lm(input.data[,1] ~ x))))
  coefficients <- cbind(intercept=scales[,1],scales[,2]*coefficients)


  best.chain <- switch(tolower(fit.metric),
                  bic = sort(fit.out[,1], index.return = TRUE)$ix,
                  tau = sort(fit.out[,1], index.return = TRUE)$ix,
                  aic = sort(fit.out[,1], index.return = TRUE)$ix
                  )

  best.coef <- coefficients[best.chain,,drop=FALSE][1,]
  fitted.values <- cbind(intercept=1,input.data[,-1]) %*% matrix(best.coef, ncol = 1)

  if (roe) {
  	roe.df <- data.frame(roe.mat)
    names(roe.df) <- c("fit.metric", "r", colnames(input.data)[-1])
  	roe.df$beta <- factor(rep(1:n.beta, times = n.gens * n.chains))
  	roe.df$gens <- factor(rep(1:n.gens, each = n.beta, times = n.chains))
  	roe.df$chain <- factor(rep(best.chain, each = n.beta * n.gens))
  }
  sim.results <- list(date = date(),
    call = match.call(),
    coefficients = coefficients[best.chain, , drop = FALSE],
    fitted.values = fitted.values,
    residuals = unlist(input.data[,1] - fitted.values),
    rank.residuals = (rank(input.data[,1]) -
        rank(fitted.values)),
    bic = c(fit.out.bic)[best.chain],
    r = c(fit.out.r)[best.chain],
    tau = c(fit.out.tau)[best.chain],
    tau.a = c(fit.out.tau.a)[best.chain],
    tau.b = c(fit.out.tau.b)[best.chain],
    tau.par = fit.out.tau.par[best.chain,],
    aic = c(fit.out.aic)[best.chain],
    metric.betas = metricbeta,
    p.vals = p.vals,
    rank = sum(best.coef != 0),
    model = data.frame(input.data),
    fit.metric = fit.metric)
  if (check.convergence) {
    sim.results$converge.fit.metric <- converge.fit.metric
    sim.results$converge.beta <- converge.beta
    sim.results$converge.r <- converge.r
    attr(sim.results, "converge.check") <- TRUE
  } else {
    attr(sim.results, "converge.check") <- FALSE
  }
  if (roe) {
    sim.results$roe <- roe.df
    attr(sim.results, "roe") <- TRUE
  } else {
    attr(sim.results, "roe") <- FALSE
  }
  if (save.results) {
    save(sim.results, file = paste(output, ".Rdata"))
  }
  return(sim.results)
}

##### Package functions #####

gemm <- function(x, ...) UseMethod("gemm")

gemm.default <- function(x, k.pen, parallel = FALSE, n.chains = 4, fit.metric = "bic",...) {
  if(!(parallel)) {
    est <- gemmEst(input.data = x, k.pen = k.pen, n.chains = n.chains, fit.metric=fit.metric,...)
    class(est) <- "gemm"
    est
  } else {
    registerDoMC()
    gemm.list <- foreach(i = 1:n.chains) %dopar% {
      gemmEst(input.data = x, k.pen = k.pen, n.chains = 1, fit.metric=fit.metric,...)
    }
    # Make everything gemm objects
    lapply(gemm.list, "class<-", "gemm")
    # Order Chains
    switch(tolower(fit.metric),
        bic = fit.name <- "bic",
        aic = fit.name <- "aic",
        tau = fit.name <- "tau"
        )
    if (fit.name == "tau") {
      chain.order <- order(sapply(gemm.list[1:n.chains], function(x) get(fit.name,x)), decreasing = TRUE)
    } else {
      chain.order <- order(sapply(gemm.list[1:n.chains], function(x) get(fit.name,x)))
    }
    gemm.ordered <- gemm.list[rank(chain.order)]
    # Select Best Chain
    best.chain <- gemm.ordered[[1]]
    best.chain$coefficients <-  do.call(rbind,lapply(gemm.ordered[1:n.chains],coefficients))
    best.chain$bic <-  sapply(gemm.ordered,function(x) x$bic)
    best.chain$r <-  sapply(gemm.ordered,function(x) x$r)
    best.chain$tau <-  sapply(gemm.ordered,function(x) x$tau)
    best.chain$tau.a <-  sapply(gemm.ordered,function(x) x$tau.a)
    best.chain$tau.b <-  sapply(gemm.ordered,function(x) x$tau.b)
    best.chain$tau.par <-  t(sapply(gemm.ordered,function(x) x$tau.par))
    best.chain$aic <-  sapply(gemm.ordered,function(x) x$aic)
    if (attr(gemm.ordered[[1]], "converge.check")) {
      best.chain$converge.fit.metric <- sapply(gemm.ordered,function(x) x$converge.fit.metric)
      best.chain$converge.beta <- sapply(gemm.ordered,function(x) x$converge.beta)
      best.chain$converge.r <- sapply(gemm.ordered,function(x) x$converge.r)
    }
    return(best.chain)
  }
}

print.gemm <- function(x, ...) {
 # Select correct fit value 
  switch(tolower(x$fit.metric),
         bic = fit <- x$bic,
         aic = fit <- x$aic,
         tau = fit <- x$tau
         )
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)

  if(x$fit.metric=="tau") {

    x$fit.metric <- ifelse(as.list(x$call["isTauB"])[[1]],"tau-b","tau-a")

  }

  cat("\n",x$fit.metric,"\n",sep="")
  print(fit)
}

gemm.formula <- function(formula, data=list(),...) {
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)[,-1]
  y <- matrix(model.response(mf), ncol = 1, dimnames = list(NULL, names(mf)[1]))
  #main effect variables
  me <- attributes(attributes(mf)$terms)$term.labels[attributes(attributes(mf)$terms)$order==1]
  mm <- model.matrix(attr(mf, "terms"), data=mf)
  fmla <- as.formula(paste(names(mf)[1], " ~ ", paste(me, collapse= "*")))

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
  names.betas.all <- names
  names.betas.in.model <- attr(terms(mf),"term.labels")
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

plot.gemm <- function(x, ...) {
  par(mfrow = c(1,3))
  if (attr(x, "converge.check")) {
    par(mfrow = c(2,2))
    convergencePlot(x$converge.fit.metric,x$fit.metric)
  }
  plot(rank(fitted.values(x)), rank(x$model[1]),
    main = "Ordinal model predictions", xlab = "Rank predictions",
    ylab = "Rank criterion")
  plot(fitted.values(x), unlist(x$model[1]), main = "Metric model predictions",
    xlab = "Predictions", ylab = "Criterion")
  plot(order(x$model[1]), x$rank.residuals[order(x$model[1])],
    main = "Rank disparity by criterion rank",
    xlab = "Ordered criterion", ylab = "Rank disparity")
    par(mfrow=c(1,1))
}

convergencePlot <- function(beta, fit.metric, ...) {
  y.lab <- switch(tolower(fit.metric),
                  bic = "BIC",
                  tau = "1 - tau",
                  aic = "AIC"
                  )

  chains <- ncol(beta)
  max.rep <- nrow(beta)
  xrange <- c(1, max.rep) 
  yrange <- c(min(beta), max(beta))
  plot(xrange, yrange, type="n", xlab = "rep #", ylab= y.lab)
  colors <- rainbow(chains) 
  linetype <- c(1:chains) 
  plotchar <- seq(1:chains)
  for (i in 1:chains) {
    lines(1:max.rep, beta[,i], type="b", lwd=1.5,
      lty=linetype[i], col=colors[i], pch=plotchar[i]) 
  }  
  title(paste("Convergence of ",y.lab,sep=""))
  legend(xrange[1], yrange[2], 1:chains, cex=0.8, col=colors, pch=plotchar,
    lty=linetype, title="Chains")
}

predict.gemm <- function(object, newdata = NULL, ..., tie.struct = FALSE) {
 if (is.null(newdata)) {
    out <- fitted(object)
  } else {
    tt <- terms(object$formula)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    beta <- object$coefficients[1,]
    out <- X %*% beta
  }
  out <- as.vector(out)
  attr(out, "tie.struct") <- data.frame(tauTest(unlist(model.frame(object)[1], use.names = FALSE), fitted(object), length(fitted(object))))
  return(out)
}

summary.gemm <- function(object, ...) {
  y <- object
  y$logLik <- logLik(y)
  y$AIC <- AIC(y)
  y$r.squared <- y$r[1]^2
  y$adj.r.squared <- 1 - (1 - y$r.squared) * (nobs(y) - 1)/(nobs(y) - sum(y$best.coef != 0) - 1)
  class(y) <- "summary.gemm"
  return(y)
}

print.summary.gemm <- function(x, ...) {
  print.gemm(x, ...)
}

logLik.gemm <- function(object, ...) {
  res <- object$residuals
  p <- sum(coefficients(object)[1,] != 0)
  N <- length(res)
  w <- rep.int(1, N)
  N0 <- N
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}

deviance.gemm <- function(object, ...) {
  return(sum(weighted.residuals(object)^2))
}

nobs.gemm <- function(object, ...) {
  return(nrow(residuals(object)))
}