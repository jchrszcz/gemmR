#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector genAlg(NumericMatrix metricbeta, int nbeta,
  int nsuperelites, int p, int reps,
  NumericMatrix bestmodels, bool seedmetric) {
  
  if (!seedmetric) {
    metricbeta = runif(p);
  }
  if (reps == 1) {
    NumericMatrix betas(p, nbeta);
    NumericVector temprand(p);
    NumericVector tempnorm(p);
    double scaling = sqrt(.1);
    double scaling2 = sqrt(.01);
    double scaling3 = sqrt(.05);
    betas(_,0) = mectricbeta;
    for (int i = 1; i < n.beta; i++) {
      temprand = runif(p); # outside loop
      tempnorm = rnorm(p); # outside loop
      if (i > 0 && i < 999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta
      }
      if (i >= 999 && i < 2999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * scaling;
      }
      if (i >= 2999 && i < 5999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * scaling2;
      }
      if (i >= 5999) {
        betas(_,i) = ifelse(temprand < .25, 1, ifelse(temprand > .75, -1,
          temprand));
        betas(_,i) = betas(_,i) * metricbeta + tempnorm * scaling3;
      }
    }
    for (int i = 0; i < nbeta; i++) {
      if (sum(betas(_,i)) == 0) {
        betas(_,i) = rnorm(p);
        temprand = runif(p); # move outside loop
        betas(_,i) = ifelse(temprand > .5, betas(_,i) * scaling, betas(_,i) * scaling * -1);
      }
    }
  }
  if (reps > 1) {
    int x = bestmodels.ncol();
    NumericMatrix superelites(p, nsuperelites);
    superelites = bestmodels(_,nsuperelites);
    NumericMatrix betasa = bestmodels(Range(1,p-1),_);
    NumericMatrix betasb = bestmodels(Range(1,p-1),_);
    int parent1;
    int parent2;
    int tempsize = floor(nbeta/2 + .5) - floor(nsuperelites/2 + .5);
    NumericVector temprand = runif(tempsize);
    NumericVector k = floor(runif(tempsize,1,5) + .5);
    NumericMatrix newx1(p,tempsize);
    NumericMatrix newx2(p,tempsize);
    for (int i = 0; i < tempsize; i++) {
      parent1 = floor(nsuperelites * runif(1) + .5);
      do { (parent2 = floor(nsuperelites * runif(1) + .5))
        } while (parent1 == parent2);
      if (tempand(i) < .85) {
        if (k == p) {
          newX1(_,i) = betasa(_,parent1);
          newX2(_,i) = betasb(_,parent1);
        }
        if (k < p) {
          newX1(Range(0,(k-1)),i) = betasa(Range(0,(k-1)),parent1);
          newX1(Range(k,(p-1)),i) = betasb(Range(k,(p-1)),parent2);
          newX2(Range(0,(k-1)),i) = betasb(Range(0,(k-1)),parent1);
          newX2(Range(k,(p-1)),i) = betasa(Range(k,(p-1)),parent2);
        }
      }
      if (temprand(i) >= .85) {
        newX1 = betasa(_,parent1);
        newX2 = betasb(_,parent1);
      }
    }
    
/*
      temp.rand <- matrix(runif(n.beta*p), ncol = (p*2))
    temp.rand.2 <- matrix(runif(n.beta*p), ncol = (p*2))
    for (i in 1:p) {
      new.X1[,i] <- ifelse(temp.rand[,i] < .01, temp.rand[,(i + p)], new.X1[,i])
      new.X2[,i] <- ifelse(temp.rand.2[,i] < .01,
                           temp.rand.2[,(i + p)], new.X2[,i])
    }
    super.elites <- super.elites[,-1]
    betas <- rbind(as.matrix(super.elites), new.X1, new.X2)
*/    
  }
  return y;
}

/*
#if (reps > 1) {
#    size <- dim(bestmodels)
#    elites <- bestmodels
#    sorted.elites <- elites[order(elites[,1]),]
#    super.elites <- sorted.elites[1:n.super.elites,]
#    temp.betas.a <- as.matrix(sorted.elites[,2:size[2]])
#    temp.betas.b <- as.matrix(sorted.elites[,2:size[2]])
#    parent.1 <- round(1 + (size[1] - 1) * runif(1))
#    parent.2 <- 0
#    while (parent.1 == parent.2 | parent.2 == 0) {
#      parent.2 <- round(1 + (size[1] - 1) * runif(1))
#    }
#    temp.rand <- runif(n.beta/2)
#    new.X1 <- matrix(rep(0, times = (n.beta/2 * p)), ncol = p)
#    new.X2 <- new.X1
#    for (i in 1:(n.beta/2)) {
#      parent.1 <- round(1 + (size[1] - 1) * runif(1))
#      parent.2 <- 0
#      while (parent.1 == parent.2 | parent.2 == 0) {
#        parent.2 <- round(1 + (size[1] - 1) * runif(1))
#      }
#      if (temp.rand[i] < .85) {
#        k <- round(1 + ((size[2] - 1) * runif(1)))
#        if (k == p) {
#          new.X1[i,] <- as.numeric(c(temp.betas.a[parent.1,]))
#          new.X2[i,] <- as.numeric(c(temp.betas.b[parent.1,]))
#        }
#        if (k < p) {
          new.X1[i,] <- as.numeric(c(temp.betas.a[parent.1,1:k], temp.betas.b[parent.2, ((k+1):p)]))
          new.X2[i,] <- as.numeric(c(temp.betas.b[parent.1,1:k], temp.betas.a[parent.2, ((k+1):p)]))
#        }
#      }
#      if (temp.rand[i] >= .85) {
        new.X1[i,] <- as.numeric(temp.betas.a[parent.1,])
        new.X2[i,] <- as.numeric(temp.betas.b[parent.2,])
#      }
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
  NumericMatrix randTest() {
    NumericMatrix x(5);
    NumericVector a = rnorm(3);
    x(Range(0,2),1) = a;
    return x;
  }
')
/*