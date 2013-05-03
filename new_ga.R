##### New Genetic Algorithm #####

require(Rcpp)

## REMEMBER WE'RE FILLING BY COLUMNS ##
# need to transpose bestmodels on R side
##### Begin Rcpp #####
// [[Rcpp::export]]
NumericVector genAlg(NumericMatrix metricbeta, double nbeta,
  double nsuperelites, double p, double reps,
  NumericMatrix bestmodels, bool seedmetric) {
  
  if (!seedmetric) {
    metricbeta = runif(p);
  }
  if (reps == 1) {
    NumericMatrix betas(p, nbeta);
    double scaling = sqrt(.1);
    betas(_,0) = mectricbeta;
    for (int i = 1; i < n.beta; i++) {
      temprand = runif(p);
      tempnorm = rnorm(p);
      if (i > 0 && i < 999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta
      }
      if (i >= 999 && i < 2999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * sqrt(.1);
      }
    if (i >= 2999 && i < 5999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * sqrt(.01);
      }
      if (i >= 5999) {
        betas(_,i) = ifelse(temprand < .25, 1, ifelse(temprand > .75, -1,
          temprand));
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * sqrt(.05);
      }
    }
    for (int i = 0; i < nbeta; i++) {
      if (sum(betas(_,i)) == 0) {
        betas(_,i) = rnorm(p);
        temprand = runif(p);
        betas(_,i) = ifelse(temprand > .5, betas(_,i) * scaling, betas(_,i) * scaling * -1);
      }
    }
  }
  if (reps > 1) {
    int x;
    x = bestmodels.ncol();
    NumericMatrix superelites(p, nsuperelites);
    superelites = bestmodels(_,nsuperelites);
    NumericMatrix betasa(p, x);
    NumericMatrix betasb(p, x);
    betasa = bestmodels(Range(1,p-1),_);
    betasb = bestmodels(Range(1,p-1),_);
    
    
  
  }
}

#if (reps > 1) {
#    size <- dim(bestmodels)
#    elites <- bestmodels
#    sorted.elites <- elites[order(elites[,1]),]
#    super.elites <- sorted.elites[1:n.super.elites,]
#    temp.betas.a <- as.matrix(sorted.elites[,2:size[2]])
#    temp.betas.b <- as.matrix(sorted.elites[,2:size[2]])
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

##### End Rcpp #####
  
#### Test functions

cppFunction('
  NumericVector randTest() {
    NumericVector x(1);
    x = floor(runif(1, 1, 4) + .5);
    return x;
  }
')