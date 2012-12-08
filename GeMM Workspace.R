setwd("/Users/mikedougherty/Chrabaszcz/Dropbox/Learn R!/GeMM")
setwd("/home/jcz/Dropbox/Learn R!/GeMM")
setwd("J:/my dropbox/Learn R!/GeMM")

data <- read.csv("Rgemmtest.csv")[-1]
data <- read.csv("cultureofhonor_orig92.csv")

source("gemmquick dev.R")
source("gemmquick.R")
source("genopt.R")
test <- gemmModel(data, "Rgemmoutput.Rdata", 2000, .5, 0, "orig", 0, 10, 10)

test <- gemmModel(data, "Rgemmoutput.Rdata", 500, .5, 0, "orig", 0, 3, 3)