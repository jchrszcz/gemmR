/*
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
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix genAlg(NumericMatrix metricbeta, int nbeta,
  int nsuperelites, int p, int gens,
  NumericMatrix bestmodels, bool seedmetric) {
  if (!seedmetric) {
    NumericVector tempmetric = runif(p);
    for (int i = 0; i < p; i++) {
      metricbeta(i,0)= tempmetric(i);
    }
  }
  
  NumericMatrix betas(p, nbeta);
 
  if (gens == 1) {
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
  
  if (gens > 1) {
    // assign beta vectors from bestmodels to superelites
    NumericMatrix superelites(p, nsuperelites);
    superelites = bestmodels(_,Range(0,(nsuperelites-1)));
    
    // why isn't it bestmodels(Range(1,p),_)?
    NumericMatrix betasa = bestmodels(Range(1,1 + (p-1)),_);
    NumericMatrix betasb = bestmodels(Range(1,1 + (p-1)),_);
    
    //
    NumericVector parent1(1);
    NumericVector parent2(1);
    
    // tempsize = half of # betas - half of # superelites
    int tempsize = floor(nbeta/2 + .5) - floor(nsuperelites/2 + .5);
    NumericVector temprand = runif(tempsize);
    
    // k = tempsize vector of random integers {1,2,3,4,5}
    NumericVector k = floor(runif(tempsize,1,p) + .5);
    
    //
    NumericMatrix newx1(p,tempsize);
    NumericMatrix newx2(p,tempsize);
    for (int i = 0; i < tempsize; i++) {

      // pick a random integer {1,...,nsuperelites}
      parent1 = floor(nsuperelites * runif(1) + .5);
      
      // pick a different random integer {1,...,nsuperelites}
      do { (parent2 = floor(nsuperelites * runif(1) + .5));
        } while (parent1(0) == parent2(0));

      //
      // 85% of the time
      if (temprand(i) < .85) {
        // copy parameter column 
        // if number of predictors is equal to random selection from {1,2,3,4,5} ???
        // betasa == betasb, so why use both?
        if (k(i) == p) {
          newx1(_,i) = betasa(_,parent1(0));
          newx2(_,i) = betasb(_,parent1(0));
        }
        // 
        if (k(i) < p) {
          for (int j = 0; j < p; j++) {
            if (j < k(i)) {
              newx1(j,i) = betasa(j,parent1(0));
              newx2(j,i) = betasb(j,parent1(0));
            }
            if (j >= k(i)) {
              newx1(j,i) = betasa(j,parent2(0));
              newx2(j,i) = betasb(j,parent2(0));            
            }
          }
        }
        // Can't k(i) be > p ??
      }
      
      //
      // 15% of the time
      if (temprand(i) >= .85) {
        newx1(_,i) = betasa(_,parent1(0));
        newx2(_,i) = betasb(_,parent1(0));
      }
    }
    
    //
    NumericVector temprand3 = runif(p*tempsize);
    NumericVector temprand2 = runif(p*tempsize);
    
    // nice trick to treat vector as matrix by the way.
    // 1% chance to change individual parameter value by runif(1)
    for (int i = 0; i < tempsize; i++) {
      for (int j = 0; j < p; j++) {
        if (temprand3(j + (p*i)) < .01) {
          NumericVector tempnorm = newx1(j,i)*rnorm(1,1,.25);
          newx1(j,i) = tempnorm(0);
        }
        if (temprand2(j + (p*i)) < .01) {
          NumericVector tempnorm = newx2(j,i)*rnorm(1,1,.25);
          newx2(j,i) = tempnorm(0);
        }
      }
    }
    
    // Where do you use newx2()?
    
    superelites = superelites(Range(1,(1 + (p-1))),_);
    for (int i = 0; i < nbeta; i++) {
      if (i < nsuperelites) {
        betas(_,i) = superelites(_,i);
      } else if (i < (nsuperelites + tempsize)) {
        betas(_,i) = newx1(_,(i-nsuperelites));
      } else {
        betas(_,i) = newx2(_,(i-nsuperelites-tempsize));
      }
    }
  }
  NumericMatrix allbetas = betas;
  return allbetas;
}
