## ----setup,message=FALSE,tidy=TRUE---------------------------------------
library(gemmR)
data(culture)
mod <- gemm(murder.rate ~ pasture + gini + gnp, data = culture, n.chains = 3, n.gens = 10, n.beta = 200, check.convergence = TRUE)

## ------------------------------------------------------------------------
summary(mod)

## ----plotting,tidy=TRUE--------------------------------------------------
plot(mod)

## ------------------------------------------------------------------------
yhat <- predict(mod, tie.struct = TRUE)
head(yhat)
attr(yhat, "tie.struct")

## ------------------------------------------------------------------------
logLik(mod)
AIC(mod)
BIC(mod)

