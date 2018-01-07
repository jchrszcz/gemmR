gemmR
=====

General Monotone Model in R. This is an R version of the non-parametric regression algorithm proposed by Dougherty & Thomas (2012). The algorithm works as follows:

1. Generate random weights with a genetic algorithm
2. Test the fit between the sum of the weighted predictors and an outcome using Bayesian Information Criterion based on tau-to-r transformation.
3. Select top 1/16th of candidate beta vectors, mutate and generate new weights, repeat for some number of generations.
4. Repeat entire process some number of times to check for convergence.

The resulting `gemm` object is an S3 class object.
We're working on increasing usability and writing a proper vignette.
Information on the basis for GeMM can be found in the [original paper](http://www.bsos.umd.edu/psyc/dougherty/pdf%20articles/DoughertyThomas2012Rev.pdf) or the [followup](http://onlinelibrary.wiley.com/doi/10.1111/bmsp.12090/full) [papers](http://damlab.umd.edu/pdf%20articles/2015%20Sociological%20Methodology-2015-Dougherty-223-71.pdf).

Installation
-----

We recommend installing with `devtools`:


```r
library(devtools)
install_github("jchrszcz/gemmR")
```

`gemmR` requires `Rcpp`, which also means you'll need a C++ compiler. Clear and reliably-updated directions for installing and troubleshooting those things are maintained at [the STAN project page](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#prerequisites).

Parallel Usage
------

To make `gemmR` more consistent with other R packages, we've removed the `parallel` argument.
Instead, we've added a helper function to take the output of a parallelization process (which might differ based on your preference and OS) and produce a `gemm` object.
Here's an example, if you're using windows you might need to replace `doMC` with `doSnow` and `%dopar%` with `%do%`:


```r
library(doMC)
registerDoMC()
fit <- foreach (i = 1:3) %dopar% {
      gemm(mpg ~ disp + cyl, data = mtcars, n.chains = 1, n.gens = 3, n.beta = 200)
    }
gemm.model <- list2gemm(fit)
```

New Version 1.3.03 (1-07-18)
------

* gave up on CRAN
* added C++11 flag to allow for compilation on Windows


Licensing
-----

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
