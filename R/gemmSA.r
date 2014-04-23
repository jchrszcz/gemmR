##### gemmSA TODO
# - OCLO

require(GenSA)
require(gemmR)
require(Rcpp)
require(bestglm)

fitFun <- function(x, dat) {
  y <- dat[,1]
  X <- dat[,-1, drop = FALSE]
  return(1 - tauTest(y, X %*% matrix(x, ncol = 1), length(y))$"tau-b")
}

fitGemm <- function(dat, k, minimum = -1, maximum = 1, tol = 1e-16, global.min = 0, ...) {
  lower <- rep(minimum, k)
  upper <- rep(maximum, k)
  out <- GenSA(lower = lower, upper = upper, fn = fitFun, control=list(threshold.stop=global.min+tol), dat = dat)
  # joe's scaling trick
  scales <- coef(lm(dat[,1] ~ matrix(dat[,-1], ncol = k) %*% matrix(out$par, ncol = 1)))
  coefficients <- c(intercept=scales[1],scales[-1]*out$par)
  return(list(coefficients = coefficients, fit = 1 - out$value, k = k))
}

bicTau <- function(tau, p, k, n) {
  rp <- sin(pi/2 * tau * (n - p - 1)/n)
  return(n * log(1 - rp^2) + k * log(n))
}

subFun <- function(x, y, X, p) {
    vars <- as.logical(rev(to.binary(x, p)))
    k <- sum(vars)
    Xi <- X[, vars, drop = FALSE]
    tmp <- fitGemm(cbind(y, Xi), k)
    vars[vars != 0] <- tmp$coefficients[-1]
    coefs <- c(tmp$coefficients[1], vars)
    return(c(coefs, k, tmp$fit, cor(y, cbind(1,X) %*% coefs), bicTau(tmp$fit, p, k, nrow(y))))
}

gemmFun <- function(dat, parallel = FALSE, ...) {
  p <- ncol(dat) - 1
  M <- 2^p - 1
  y <- dat[, 1, drop = FALSE]
  X <- dat[, -1, drop = FALSE]
  if (parallel) {
      require(parallel)
      all.models <- mclapply(1:M, function(x) subFun(x, y, X, p))
    } else {
    all.models <- lapply(1:M, function(x) subFun(x, y, X, p))
  }
  all.models <- simplify2array(all.models)
  all.models <- data.frame(t(cbind(all.models, c(mean(y), rep(0, p), rep(0, 4)))))
  names(all.models) <- c("(intercept)", colnames(dat)[-1], "k", "tau", "r", "bicT")
  return(all.models)
}

##### test

x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
x4 <- rnorm(100)
y <- x1 + x2 + x4
dat <- cbind(y, x1, x2, x3)

system.time(temp1 <- gemmFun(dat))
system.time(temp2 <- gemmFun2(dat))
system.time(temp3 <- gemm(y ~ x1 + x2 + x3, data = data.frame(dat), fit.metric = "tau", n.chains = 1))