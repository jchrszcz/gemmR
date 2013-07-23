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

//' genAlg
//' 
//' Generatic algorithm for gemmR.
//' @param metricbeta Weights derived using multiple regression. Overwritten if \code{seedmetric} is FALSE.
//' @param nbeta Number of candidate weight vectors in each generation.
//' @param nsuperelites Number of candidate weight vectors to involve in permutation for reps > 1.
//' @param p Number of potential predictors.
//' @param reps Generation number. For reps == 1, entirely new weights are generated. When reps > 1, bestmodels are used to generate permutations.
//' @param bestmodels Matrix of best candidate weight vectors from previous generation.
//' @param seedmetric If TRUE, multiple regression weights are used to seed the genetic algorithm. Otherwise, random weights are used.
//' @return y Matrix of candidate weights with number of rows equal to p.
//' 
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix genAlg(NumericMatrix metricbeta, int nbeta, int nsuperelites, int p, int reps, NumericMatrix bestmodels, bool seedmetric) {
  if (!seedmetric) {
    NumericVector tempmetric = runif(p);
    for (int i = 0; i < p; i++) {
      metricbeta(i,0)= tempmetric(i);
    }
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
    NumericMatrix betasa = bestmodels(Range(1,1 + (p-1)),_);
    NumericMatrix betasb = bestmodels(Range(1,1 + (p-1)),_);
    NumericVector parent1(1);
    NumericVector parent2(1);
    int tempsize = floor(nbeta/2 + .5) - floor(nsuperelites/2 + .5);
    NumericVector temprand = runif(tempsize);
    NumericVector k = floor(runif(tempsize,1,5) + .5);
    NumericMatrix newx1(p,tempsize);
    NumericMatrix newx2(p,tempsize);
    for (int i = 0; i < tempsize; i++) {
      parent1 = floor(nsuperelites * runif(1) + .5);
      do { (parent2 = floor(nsuperelites * runif(1) + .5));
        } while (parent1(0) == parent2(0));
      if (temprand(i) < .85) {
        if (k(i) == p) {
          newx1(_,i) = betasa(_,parent1(0));
          newx2(_,i) = betasb(_,parent1(0));
        }
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
      }
      if (temprand(i) >= .85) {
        newx1(_,i) = betasa(_,parent1(0));
        newx2(_,i) = betasb(_,parent1(0));
      }
    }
    NumericVector temprand3 = runif(p*tempsize);
    NumericVector temprand2 = runif(p*tempsize);
    for (int i = 0; i < tempsize; i++) {
      for (int j = 0; j < p; j++) {
        if (temprand3(j + (p*i)) < .01) {
          NumericVector tempnorm = runif(1);
          newx1(j,i) = tempnorm(0);
        }
        if (temprand2(j + (p*i)) < .01) {
          NumericVector tempnorm = runif(1);
          newx2(j,i) = tempnorm(0);
        }
      }
    }
    superelites = superelites(Range(1,(1 + (p-1))),_);
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
