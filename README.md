gemmR
=====

General Monotone Model in R. This is an R version of the non-parametric regression algorithm proposed by Dougherty & Thomas (2012). The algorithm works as follows:

1. Generate random weights with a genetic algorithm
2. Test the fit between the sum of the weighted predictors and an outcome using Bayesian Information Criterion based on tau-to-r transformation.
3. Select top 1/16th of candidate beta vectors, mutate and generate new weights, repeat for some number of generations.
4. Repeat entire process some number of times to check for convergence.

The resulting `gemm` object is currently an S3 class object. We're working on reducing runtime, (mostly through parallelization at this point), and writing a proper helpfile. Information on the basis for GeMM can be found in the [original paper](http://www.bsos.umd.edu/psyc/dougherty/pdf%20articles/DoughertyThomas2012Rev.pdf).

**Bug Notes**

*  On 07/22/2013, we implemented a bugfix that allows the code to be compiled on Windows and Mac OS X. You can now install the gemmR package by downloading the `.tar.gz` file in this repository, including it in your working directory, and running either:

    `install.packages("gemmR_1.0.tar.gz", repos = NULL, type = "source")`
    
  or, if you have `devtools` installed:
    
    `install_github("gemmR", "jchrszcz", subdir = "gemmR")`

*  Predictor names that are entirely contained within other predictor names can cause problems with the `gemm.formula` function.

* `roe` argument in `gemmEst` is not currently functioning.

**Licensing Note**

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
