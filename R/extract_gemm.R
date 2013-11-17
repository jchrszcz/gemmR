library(texreg)

# add recommend flag for texreg

#standard (bestcoef, log(L)?, aic, nobs, r2)
#gemm (coef, bic, r2, tau, nobs)
#cross (estcoef, crosscoef, r tau bic)
#report (bestcoef, bic, tau, R2, n)

#remove best chain
#coefficients
#est bic r tau
#cross bic r tau
#standard model fits (LL? and AIC)

extract.gemm <- function(model,include.rsquared = TRUE, include.tau = TRUE, include.bic = TRUE, include.adjrs = FALSE, include.nobs = FALSE, include.aic = FALSE, include.lik = FALSE, include.dev = FALSE, ...) {
  s <- summary(model, ...)
  names <- colnames(coef(s))
  co <- s$best.coef
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared) {
    rs <- s$r.squared
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs) {
    adj <- s$adj.r.squared
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs) {
    n <- nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.aic) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.lik) {
    lik <- logLik(model)
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dev) {
    dev <- deviance(model)
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.tau) {
    tau <- model$est.tau[1]
    gof <- c(gof, tau)
    gof.names <- c(gof.names, "$\tau$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  tr <- createTexreg(
    coef.names = names,
    coef = co,
    se = rep(0, times = length(co)),
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
  
  return(tr)
}

setOldClass("gemm", test = FALSE, className("gemm", "gemmR"))

setMethod("extract", signature = className("gemm", package = "gemmR"), definition = extract.gemm)

# example code
#newt <- extract(d)
#{ sink("/dev/null"); string <- texreg(newt, return.string = TRUE); sink(); }
#list.o.lines <- strsplit(string,"\\n",useBytes=TRUE)[[1]]
#cat(list.o.lines[!grepl("0)$",list.o.lines,fixed=TRUE)],sep="\n")
