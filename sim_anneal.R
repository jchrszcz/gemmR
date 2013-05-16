##### Simulated Annealing #####

require(GenSA)

gemmFit <- function(n, betas, data, p, k.cor, pearson) {
  if (sum(betas == 0) == p) {
    tau <- 0
    r <- 0
  }
  if (sum(betas == 0) != p) {
    tau <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)),
               method = "kendall")
    if (pearson) {
      r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
    }
  }
  if (tau < 0) {
    betas <- betas * -1
    tau <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)),
               method = "kendall")
    if (pearson) {
      r <- cor(c(data[,1]), c(.rowSums(t(betas * t(data[,-1])), n, p)))
    }
  }
  k <- sum(betas != 0)
  knp <- sin(pi/2*tau*((n-k-1)/n))
  bic <- n * log(1 - knp ^ 2) + k * log(n)
  y <- list(bic = bic)
  if (pearson) {
    y$r <- r
    y$tau <- tau
  }
  return(y)
}

##### Example #####

# Try Rastrgin function (The objective function value for global minimum
# is 0 with all components of par are 0.)
Rastrigin <- function(x) {
sum(x^2 - 10 * cos(2 * pi * x)) + 10 * length(x)
}
set.seed(1234) # The user can use any seed.
dimension <- 30
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = Rastrigin,
  control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par","counts")]
# GenSA will stop after running for about 2 seconds
# Note: The time for solving this problem by GenSA may vary
# depending on the computer used.
set.seed(1234) # The user can use any seed.
dimension <- 30
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = Rastrigin,
control=list(max.time=2))
out[c("value","par","counts")]

##### Testing #####

intercept <- rep(1, times = 100)
x1 <- rnorm(100)
x2 <- rnorm(100)
y <- as.vector(engle$ravens)
x <- as.matrix(engle[,3:ncol(engle)])

exampleFun <- function(input,x,y) {
  if (sum(input == 0) == length(input)) {
    bic <- 0
  }
  if (sum(input == 0) != length(input)) {
    prods <- x
    k <- sum(input != 0)
    for (i in 1:ncol(x)) {
      prods[,i] <- x[,i] * input[i]
    }
    pred <- rowSums(prods)
    r <- sin(.5*pi*cor(y,pred, method = "kendall"))
    bic <- length(y) * log(1 - r ^ 2) +  k*log(length(y))
  }
  return(bic)
}


dimension <- 9
lower <- rep(0, times = 9)
upper <- rep(1, times = 9)
par <- as.vector(coefficients(lm(y ~ x))[-1])
control = list(max.time = 100)
system.time(out <- GenSA(par = par, lower = lower, upper = upper, fn = exampleFun, x = x, y = y))
out[c("value","par","counts")]