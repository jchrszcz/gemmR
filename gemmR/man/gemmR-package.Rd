\name{gemmR-package}
\alias{gemmR-package}
\alias{gemmR}
\docType{package}
\title{
  General Monotone Model in R (gemmR)
}
\description{
  The \code{gemmR} package provides functions for fitting and analyzing General
  Monotone Models (GeMM). The GeMM is a statistical algorithm for predicting
  rank orders from a set of \code{k} predictors. As shown in Dougherty and
  Thomas (2012, Psychological Review), GeMM is unaffected by any monotonic
  transformation of the criterion variable, unaffected by non-linear monotone
  relationships, and shows better power and predictive accuracy than
  Least-Squares regression procedures in a variety of contexts. The primary
  \code{gemmR} function is \code{\link{gemm}}.
}
\details{
  \tabular{ll}{
    Package: \tab gemmR\cr
    Type: \tab Package\cr
    Version: \tab 1.3.01\cr
    Date: \tab 2014-09-08\cr
    License: \tab GPL-3\cr
  }
}
\author{
  Jeffrey S. Chrabaszcz\cr
  Joe Tidwell\cr
  Michael R. Dougherty\cr
}
\section{license}{
  This file is part of gemmR. gemmR is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  gemmR is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  gemmR.  If not, see \url{<http://www.gnu.org/licenses/>}.
}
\references{
  Dougherty, M. R., & Thomas, R. P. (2012). Robust decision making in a
  nonlinear world. \emph{Psychological Review}, 119(2), 321.
  \url{http://www.bsos.umd.edu/psyc/dougherty/pdf\%20articles/
  DoughertyThomas2012Rev.pdf}
}
\keyword{
  package
}
