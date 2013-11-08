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
along with gemmR.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
const float PIE = 4.0*atan(1.0);

/* Comparator function for qsortDoubles.*/
static int qsortDoubleComp(const void* lhs, const void* rhs) {
    double lhsNum, rhsNum;
    lhsNum = *((double*) lhs);
    rhsNum = *((double*) rhs);

    if(lhsNum < rhsNum) {
        return -1;
    } else if(lhsNum > rhsNum) {
        return 1;
    } else {
        return 0;
    }
}


/* Wrapper for C's qsort functionality to simplify use in this library.*/
static void qsortDoubles(double* arr, size_t len) {
    qsort(arr, len, sizeof(double), &qsortDoubleComp);
}

typedef struct {
    double first;
    double second;
} PairOfDoubles;

/* Comparator function for zipSort().*/
static int zipSortComp(const void* lhs, const void* rhs) {
    PairOfDoubles l, r;
    l = *((PairOfDoubles*) lhs);
    r = *((PairOfDoubles*) rhs);

    if(l.first < r.first) {
        return -1;
    } else if(l.first > r.first) {
        return 1;
    } else {
        return 0;
    }
}

/* Uses C's qsort functionality to sort two arrays in lockstep, ordered by
 * the first array.
 */
static void zipSort(double* arr1, double* arr2, size_t len) {
    size_t i;

    PairOfDoubles* pairs = (PairOfDoubles*) malloc(len * sizeof(PairOfDoubles));
    for(i = 0; i < len; i++) {
        pairs[i].first = arr1[i];
        pairs[i].second = arr2[i];
    }

    qsort(pairs, len, sizeof(PairOfDoubles), &zipSortComp);

    // Copy the results back.
    for(i = 0; i < len; i++) {
        arr1[i] = pairs[i].first;
        arr2[i] = pairs[i].second;
    }

    free(pairs);
}

/* Sorts in place, returns the bubble sort distance between the input array
 * and the sorted array.
 */
static uint64_t insertionSort(double* arr, size_t len) {
    size_t maxJ, i;
    uint64_t swapCount = 0;

    if(len < 2) {
        return 0;
    }

    maxJ = len - 1;
    for(i = len - 2; i < len; --i) {
        size_t j = i;
        double val = arr[i];

        for(; j < maxJ && arr[j + 1] < val; ++j) {
            arr[j] = arr[j + 1];
        }

        arr[j] = val;
        swapCount += (j - i);
    }

    return swapCount;
}

static uint64_t merge(double* from, double* to, size_t middle, size_t len) {
    size_t bufIndex, leftLen, rightLen;
    uint64_t swaps;
    double* left;
    double* right;

    bufIndex = 0;
    swaps = 0;

    left = from;
    right = from + middle;
    rightLen = len - middle;
    leftLen = middle;

    while(leftLen && rightLen) {
        if(right[0] < left[0]) {
            to[bufIndex] = right[0];
            swaps += leftLen;
            rightLen--;
            right++;
        } else {
            to[bufIndex] = left[0];
            leftLen--;
            left++;
        }
        bufIndex++;
    }

    if(leftLen) {
        memcpy(to + bufIndex, left, leftLen * sizeof(double));
    } else if(rightLen) {
        memcpy(to + bufIndex, right, rightLen * sizeof(double));
    }

    return swaps;
}

/* Sorts in place, returns the bubble sort distance between the input array
 * and the sorted array.
 */
static uint64_t mergeSort(double* x, double* buf, size_t len) {
    uint64_t swaps;
    size_t half;

    if(len < 10) {
        return insertionSort(x, len);
    }

    swaps = 0;

    if(len < 2) {
        return 0;
    }

    half = len / 2;
    swaps += mergeSort(x, buf, half);
    swaps += mergeSort(x + half, buf + half, len - half);
    swaps += merge(x, buf, half, len);

    memcpy(x, buf, len * sizeof(double));
    return swaps;
}

static uint64_t getMs(double* data, size_t len) {  /* Assumes data is sorted.*/
    uint64_t Ms = 0, tieCount = 0;
    size_t i;

    for(i = 1; i < len; i++) {
        if(data[i] == data[i-1]) {
            tieCount++;
        } else if(tieCount) {
            Ms += (tieCount * (tieCount + 1)) / 2;
            tieCount++;
            tieCount = 0;
        }
    }
    if(tieCount) {
        Ms += (tieCount * (tieCount + 1)) / 2;
        tieCount++;
    }
    return Ms;
}
/* This function calculates the Kendall Tau B on a pair of C-style "arrays".
 * Note that it will completely overwrite arr1 and sort arr2, so these need
 * to be duplicated before passing them in.
 */
double kendallNlogN(double* arr1, double* arr2, size_t len) {
    uint64_t m1 = 0, m2 = 0, tieCount, swapCount, nPair;
    int64_t s;
    size_t i;

    zipSort(arr1, arr2, len);
    nPair = (uint64_t) len * ((uint64_t) len - 1) / 2;
    s = nPair;

    tieCount = 0;
    for(i = 1; i < len; i++) {
        if(arr1[i - 1] == arr1[i]) {
            tieCount++;
        } else if(tieCount > 0) {
            qsortDoubles(arr2 + i - tieCount - 1, tieCount + 1);
            m1 += tieCount * (tieCount + 1) / 2;
            s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
            tieCount++;
            tieCount = 0;
        }
    }
    if(tieCount > 0) {
        qsortDoubles(arr2 + i - tieCount - 1, tieCount + 1);
        m1 += tieCount * (tieCount + 1) / 2;
        s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
        tieCount++;
    }

    // Using arr1 as the buffer because we're done with it.
    swapCount = mergeSort(arr2, arr1, len);

    m2 = getMs(arr2, len);
    s -= (m1 + m2) + 2 * swapCount;

    double denominator1 = nPair - m1;
    double denominator2 = nPair - m2;
    double cor = s / sqrt(denominator1) / sqrt(denominator2);
    return cor;
}

// [[Rcpp::export]]
double kt(NumericVector arr1, NumericVector arr2, int length) {
 std::vector<double> y(length);
 NumericVector xx(clone(arr1));
 NumericVector yy(clone(arr2));
  
return kendallNlogN(xx.begin(), yy.begin(), length);
}

// [[Rcpp::export]]
NumericVector fitValues (NumericVector betas, NumericMatrix data) {
  NumericVector values(data.nrow());
  for (int i = 0; i < data.nrow(); i++) {
    for (int j = 0; j < betas.size(); j++) {
      values(i) = values(i) + data(i,j+1)*betas(j);
    }  
  }
  return(values);
}

// [[Rcpp::export]]
double corRcpp(NumericVector x, const NumericVector y)
{
    size_t n = x.size();
    double ex(0), ey(0), sxx(0), sxy(0), syy(0), xt(0), yt(0);
    double tiny = 1e-20;

  for (size_t i = 0; i < n; i++) { // Find the means.
   ex += x[i];
   ey += y[i];
  }
  ex /= n;
  ey /= n;
  for (size_t i = 0; i < n; i++) { // Compute the correlation coeﬃcient.
    xt = x[i] - ex;
    yt = y[i] - ey;
    sxx += xt * xt;
    syy += yt * yt;
    sxy += xt * yt;
  }
return sxy/(sqrt(sxx*syy)+tiny);
}


// [[Rcpp::export]]
NumericVector gemmFitRcpp(int n, NumericVector betas, NumericMatrix data, int p, int kCor, bool pearson, bool correction) {
  double tau, r;

if (sum(betas == 0) == p) {
    tau = 0;
    r = 0;
  }

  if (sum(betas == 0) != p) {
    tau = kt(data(_,0), fitValues(betas,data), data.nrow());    
//      if(pearson) {
      r = corRcpp(data(_,0), fitValues(betas,data));
//    }
  }

//  if (tau < 0) {
//    betas = betas * -1;
//    tau = kt(data(_,0), fitValues(betas,data), data.nrow());
//    if(pearson) {
//       r = corRcpp(data(_,0), fitValues(betas,data));
//    }
//  }

  double knp, bic, bicr; // = sin( (PIE/2) * tau * ((n-kCor-1)/n));
  
  if(correction) {
    knp = sin(PIE/2*tau*(n-p-1)/n);
    bic = n * log(1 - pow(knp,2)) + kCor * log(n);
    bicr = n * log(1 - pow(r,2)) + kCor * log(n);
  } else {
    knp = sin(PIE/2*tau*(n-p-1)/n);
    bic = n * log(1 - pow(knp,2)) + p * log(n);
    bicr = n * log(1 - pow(r,2)) + p * log(n);    
  }
  
  return Rcpp::NumericVector::create(tau,r,bic,bicr,knp);
}



// [[Rcpp::export]]
List gemmFitRcppI(int n, NumericMatrix betas, NumericMatrix data, int p, NumericVector kCor, bool getR, CharacterVector correction) {

  NumericVector fit;
  NumericVector fitR(betas.nrow()), fitTau(betas.nrow()), fitBIC(betas.nrow()), fitBICr(betas.nrow());

  NumericMatrix data2(clone(data));
  
  bool corType; 
  
  /*std::string s_targetname = as<std::string>(targetname)
  */
  
  (as<std::string>(correction)=="knp") ? (corType=true) : (corType=false);


for (int i=0; i < betas.nrow(); i++) {
  fit = gemmFitRcpp(n, betas(i,_), data2, p, kCor(i), getR, corType);
  if(fit[1]<0) {
    for(int j =0; j<p;j++) {
      betas(i,j) = betas(i,j)*-1;
    }
  }

    fitBIC(i) = fit[2];
    fitBICr(i) = fit[3];
  if (getR) {
      fitR(i) = fit[1];
      fitTau(i) = fit[0];
    }
  }

  return Rcpp::List::create(Rcpp::Named("r") = fitR,
                            Rcpp::Named("bic") = fitBIC,
                            Rcpp::Named("tau") = fitTau,
                            Rcpp::Named("bic.r") = fitBICr);
}

NumericMatrix stl_sort_matrix(NumericMatrix x) {
   NumericMatrix y = clone(x);
   std::sort(y.begin(), y.end());

  return y;
}


// [[Rcpp::export]]
NumericMatrix sortIt(NumericMatrix x) {
  return stl_sort_matrix(x);
}