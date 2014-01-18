gemmCV <- function(x) {
  x1 <- rnorm(50)
  x2 <- rnorm(50)
  x3 <- rnorm(50)
  y <- as.vector(scale(x1 + x3 + rnorm(50)))
  df <- data.frame(y, x1, x2, x3)
  a <- gemm(y ~ x1 + x2 + x3, data = df[1:25,], n.chains = 2, n.beta = 1000, parallel = TRUE, p.est = 1)
  p.t <- c(cor(fitted.values(a), df$y[1:25], method = "kendall"), cor(as.matrix(df[26:50, -1], ncol = 3) %*% as.matrix(coefficients(a)[1,], ncol = 1), df$y[26:50], method = "kendall"), cor(fitted.values(a), df$y[1:25]), cor(as.matrix(df[26:50, -1], ncol = 3) %*% as.matrix(coefficients(a)[1,], ncol = 1), df$y[26:50]))
  a <- gemm(y ~ x1 + x2 + x3, data = df[1:25,], n.chains = 2, n.beta = 1000, parallel = FALSE, p.est = 1)
  p.f <- c(cor(fitted.values(a), df$y[1:25], method = "kendall"), cor(as.matrix(df[26:50, -1], ncol = 3) %*% as.matrix(coefficients(a)[1,], ncol = 1), df$y[26:50], method = "kendall"), cor(fitted.values(a), df$y[1:25]), cor(as.matrix(df[26:50, -1], ncol = 3) %*% as.matrix(coefficients(a)[1,], ncol = 1), df$y[26:50]))
  return(c(p.t, p.f))
}

a <- replicate(100, gemmCV(1), simplify = TRUE)
