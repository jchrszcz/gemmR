\name{culture}
\docType{data}
\alias{culture}
\title{Geographic and cultural predictors of country homicide}
\description{
  The \code{culture} data contains 92 rows and 5 columns from a study on the predictors of homicide rate by country.
}

\usage{culture}
\format{
  An object of class \code{data.frame} containing the following columns:

  \describe{
    \item{country}{A factor vector containing the country names correspond to each row of observations.}
    \item{murder.rate}{Numeric vector  giving murder rate as expected annual murders per 10,000 people.}
    \item{pasture}{Numeric vector giving the percentage of geographic area used as pasture.}
    \item{gini}{Numeric vector of Gini coefficients, a measure of income disparity.}
    \item{gnp}{Numeric vector of gross national products in US dollars.}
  }
}
\source{
  \enumerate{
    \item Dougherty, M.R., Thomas, R.P., Brown, R.P., Chrabaszcz, J.S., and Tidwell, J.W. (2014). An introduction to the general monotone model with application to two problematic datasets. Sociological Methods.
    \item Brown, R. P., & Osterman, L. L. (2012). Culture of Honor, Violence. The Oxford handbook of evolutionary perspectives on violence, homicide, and war, 218.
  }	
}