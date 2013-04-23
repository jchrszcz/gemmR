setwd("/Users/mikedougherty/Chrabaszcz/Dropbox/learnr/GeMM")
setwd("/home/jcz/Dropbox/learnr/GeMM")
setwd("/home/jchrszcz/Dropbox/learnr/GeMM")
setwd("J:/my dropbox/learnr/GeMM")

lapply(c("plyr","ggplot2","doMC","reshape"), require, character.only = TRUE)

library(plyr)
library(doMC)

data <- read.csv("Rgemmtest.csv")
data <- read.csv("cultureofhonor_orig92.csv")

source("gemmquick dev.R")
source("/R/gemmquick.R")
source("genopt.R")

system.time(test <- gemmModel(data, "gemm", 200, .5, 3, 3))

system.time(test <- gemmModel(data, "gemm", 6000, .5, 3, 3))

##### ddply development (probably scrap) #####

system.time(for (i in 1:dim(betas)[1]) {cor(c(data[,1]), c(.rowSums(t(betas[1,] * t(data[,-1])), n, p)))})
system.time(for (i in 1:dim(betas)[1]) {cor(c(data[,1]), c(.rowSums(t(betas[1,] * t(data[,-1])), n, p)), method = "kendall")})

system.time(for (i in 1:dim(betas)[1]) {
  gemm.fit.out <- gemmFit(n, betas[i,], data, p)
  fit.stats[i,] <- gemm.fit.out$bic
  fit.stats.r[i,] <- gemm.fit.out$r
})

registerDoMC()

betas <- geneticAlgorithm(metric.beta, n.beta, n.beta.vectors, n.super.elites, p, reps, bestmodels)

iterations <- 100
time.test <- matrix(rep(0, times = 4 * iterations), ncol = 4)
for (j in 1:iterations) {
  # time test for current functions
  time.test[j,1] <- as.numeric(system.time(for (i in 1:dim(betas)[1]) {cor(c(data[,1]), c(.rowSums(t(betas[1,] * t(data[,-1])), n, p)))})[1])
  time.test[j,2] <- as.numeric(system.time(for (i in 1:dim(betas)[1]) {cor(c(data[,1]), c(.rowSums(t(betas[1,] * t(data[,-1])), n, p)), method = "kendall")})[1])
  # time test for parallel adply
  time.test[j,3] <- as.numeric(system.time(a <- adply(betas, .margins = 1, function(x) cor(c(data[,1]), c(.colSums(c(as.numeric(x[1:p])) * t(data[,-1]), p, n))), .parallel = TRUE))[1])
  time.test[j,4] <- as.numeric(system.time(b <- adply(betas, .margins = 1, function(x) cor(c(data[,1]), c(.colSums(c(as.numeric(x[1:p])) * t(data[,-1]), p, n)), method = "kendall"), .parallel = TRUE))[1])
}

a <- melt(time.test)
a$package <- factor(rep(c("base","plyr"), each = 200))
a$stat <- factor(rep(c("r","tau"), each = 100, times = 2))
qplot(value, fill = package, geom = "density", data = a) + facet_wrap(~stat)