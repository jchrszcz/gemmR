gemmR
=====

General Monotone Model in R. This is an R version of the non-parametric regression algorithm proposed by Dougherty & Thomas (2012). The algorithm works as follows:

1. Generate random weights with a genetic algorithm
2. Test the fit between the sum of the weighted predictors and an outcome using Bayesian Information Criterion based on tau-to-r transformation.
3. Select top 1/16th of candidate beta vectors, mutate and generate new weights, repeat for some number of generations.
4. Repeat entire process some number of times to check for convergence.

The resulting `gemm` object is an S3 class object. We're working on increasing usability and writing a proper vignette. Information on the basis for GeMM can be found in the [original paper](http://www.bsos.umd.edu/psyc/dougherty/pdf%20articles/DoughertyThomas2012Rev.pdf).

Installation
-----

You can install `gemmR` by downloading the `.tar.gz` file in this repository, including it in your working directory, and running either:


```r
# to install from .tar.gz
install.packages("gemmR_1.2.01.tar.gz", repos = NULL, type = "source")

# to install directly from github
library(devtools)
install_github("gemmR", "jchrszcz", subdir = "gemmR", ref = "norming")
```


`gemmR` requires `Rcpp`, which also means you'll need a C++ compiler. Clear and reliably-updated directions for installing and troubleshooting those things are maintained at [the STAN project page](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites).

New Version 1.2.01 (3-22-14)
------
### Changes
* GeMM coefficients now metric scaled
* updated for Rcpp 0.11
* removed sum-to-1 norming for coefficients

TODO (12-28-13)
------
* Remove redundance for tau.a/tau.b
* Add tau.a/tau.b for cross-validation
* Enable plot.gemm even when check.convergence==FALSE

Bug Notes
-----

12-8-13

* ```cross.val.tau``` and ```cross.val.tau``` were flipped when performing cross-validation and using ```parallel=TRUE```.
* When ```n.gens``` was 2 or less, wasn't returning optimized values.
* Many misc minor fixes and optimization improvements.
11-18-13

* fixed check so NA values in `lin.mod` returns a warning and sets `seed.metric` to FALSE.
* `roe` functioning again
* `p.est` now returns cross-validation fits correctly, disables `roe` if enabled.

11-8-13

*  On 07/22/2013, we implemented a bugfix that allows the code to be compiled on Windows and Mac OS X.
*  Predictor names that are entirely contained within other predictor names can cause problems with the `gemm.formula` function.
* `roe` argument in `gemmEst` is not currently functioning.


Licensing
-----

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
