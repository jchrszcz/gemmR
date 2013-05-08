#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix genAlg(NumericVector metricbeta, int nbeta,
  int nsuperelites, int p, int reps,
  NumericMatrix bestmodels, bool seedmetric) {
  
  if (!seedmetric) {
    metricbeta = runif(p);
  }
  NumericMatrix betas(p, nbeta);
  if (reps == 1) {
    NumericVector temprand(p);
    NumericVector tempnorm(p);
    double scaling = sqrt(.1);
    double scaling2 = sqrt(.01);
    double scaling3 = sqrt(.05);
    betas(_,0) = metricbeta;
    for (int i = 1; i < nbeta; i++) {
      temprand = runif(p); // outside loop
      tempnorm = rnorm(p); // outside loop
      if (i > 0 && i < 999) {
        betas(_,i) = ifelse(temprand < .5, 1, 0);
        betas(_,i) = betas(_,i) * metricbeta;
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
        temprand = runif(p); // move outside loop
        betas(_,i) = ifelse(temprand > .5, betas(_,i) * scaling, betas(_,i) * scaling * -1);
      }
    }
  }
  if (reps > 1) {
    NumericMatrix superelites(p, nsuperelites);
    superelites = bestmodels(_,Range(0,(nsuperelites-1)));
    NumericMatrix betasa = bestmodels(Range(1,p-1),_);
    NumericMatrix betasb = bestmodels(Range(1,p-1),_);
    NumericVector parent1(1);
    NumericVector parent2;
    int tempsize = floor(nbeta/2 + .5) - floor(nsuperelites/2 + .5);
    NumericVector temprand = runif(tempsize);
    NumericVector k = floor(runif(tempsize,1,5) + .5);
    NumericMatrix newx1(p,tempsize);
    NumericMatrix newx2(p,tempsize);
    for (int i = 0; i < tempsize; i++) {
      parent1 = floor(nsuperelites * runif(1) + .5);
      do { (parent2(1) == floor(nsuperelites * runif(1) + .5));
        } while (parent1(1) == parent2(1));
      if (temprand(i) < .85) {
        if (k(i) == p) {
          newx1(_,i) = betasa(_,parent1(1));
          newx2(_,i) = betasb(_,parent1(1));
        }
        if (k(i) < p) {
          for (int j = 0; j < p; j++) {
            if (j < k(i)) {
              newx1(j,i) = betasa(j,parent1(1));
              newx2(j,i) = betasb(j,parent1(1));
            }
            if (j >= k(i)) {
              newx1(j,i) = betasa(j,parent2(1));
              newx2(j,i) = betasb(j,parent2(1));            
            }
          }
        }
      }
      if (temprand(i) >= .85) {
        newx1(_,i) = betasa(_,parent1(1));
        newx2(_,i) = betasb(_,parent1(1));
      }
    }
    superelites = superelites(Range(1,(p-1)),_);
    for (int i = 0; i < nbeta; i++) {
      if (i < nsuperelites) {
        betas(_,i) = superelites(_,i);
      } else if (i < (nsuperelites + tempsize)) {
        betas(_,i) = newx1(_,(i-nsuperelites));
      } else {
        betas(_,i) = newx1(_,(i-nsuperelites-tempsize));
      }
    }
  }
  NumericMatrix y = betas;
  return y;
}