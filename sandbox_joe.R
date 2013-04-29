

a <- rnorm(100)
b <- rnorm(100)
c <- rnorm(100)
y <- a + b + a*b + rnorm(100)
mod <- y~a+b


source('~/Dropbox/org/projects/Jeff/gemmR/gemmquick joe.R', echo=FALSE)
gemm(mod)


source('~/Dropbox/org/projects/Jeff/gemmR/gemmquick dev.R', echo=FALSE)
gemm(mod, n.beta=200, n.data.gen=1)
gemm(weight~Time*Diet, data=ChickWeight,n.beta=200, n.data.gen=1)

debug(gemm.formula)
gemm(weight~Time + Diet, data=ChickWeight,n.beta=200, n.data.gen=1)
undebug(gemm.formula)


summary(lm(mod))
gemm(mod)

data(ChickWeight)
summary(ChickWeight)




y <- ChickWeight$weight
a <- ChickWeight$Time
b <- ChickWeight$Diet

mod <- y~a*b

cbind(diag(3),contr.treatment(3))


mycars <- mtcars
mycars$cyl <- factor(mycars$cyl)
fit <- lm(  )

mod <- mpg ~wt+cyl

mf <- model.frame(mod, data=mycars, x=TRUE)
attributes(attributes(mf)$terms)$factor[-1,]



mod <- (mpg ~ wt+cyl, data=mycars, x=TRUE)



z.out <- zelig(y ~ x1 + x2 + x3 + as.factor(state), 
               data = mydata, model = "ls")

num.con <- 2
num.cat <- 1
cat1.lvl <- 3

diag(num.con+cat1.lvl-1)
contr.treatment(3)

diag(4)

all(crossprod(cT) == diag(4)) # TRUE: even orthonormal

ff <- factor(sample(letters[1:5], 25, replace = TRUE))
diag(nlevels(ff))[ff, ]
model.matrix( ~ ff - 1)

dfB  <- data.frame(x1=c(1,0,0), x2=factor(c('a','b','c')), x3=factor(c('a','b','c')) )
model.matrix(~x1*x2, data=dfB)[,-1]

expand.grid(~x1*x2, data=dfB)

mf <- model.frame(~x1*x2, data=dfB)

attributes(attributes(mf)$terms)$factor[-1,]


summary(lm(weight~Time*Diet, data=ChickWeight))
gemm(weight~Time*Diet, data=ChickWeight)

library(help="datasets")

#### HERE ####

a <- rnorm(100)
b <- rnorm(100)

b <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))
c <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))
d <- rnorm(100)
y <- a + b + a*b + rnorm(100)
b <- factor(b)
c <- factor(c)
mod <- y~a*b
mm <- model.matrix(mod)
mf <- model.frame(mod)

summary(lm(mod))
gemm(mod)


attributes(mm)
attributes(mf)

#main effect variables
me <- attributes(attributes(mf)$terms)$term.labels[attributes(attributes(mf)$terms)$order==1]

#full combination
count <- 0
for (var in me) {
  if(is.factor(mf[,var])) {
    print(length(levels(mf[,var]))-1)
    count <- count + length(levels(mf[,var]))-1
  print(count)
    } else {
    count = count +1
  }
}
lst <- list(NULL)
for(i in 1:count) {
  lst[[i]] <- c(0,1)  
}
lst <- t(expand.grid(lst))[,-1]
lst <- lst[,order(colSums(lst))]
rownames(lst) <- attributes(mm)$dimnames[[2]][-1][!grepl(":",attributes(mm)$dimnames[[2]][-1])]
lst <- data.frame(lst)
lst$assign <- attributes(mm)$assign[-1][1:count]
#lst[,assign:=attributes(mm)$assign[-1][1:count]]

tmp <- rep(TRUE, times=dim(lst)[2])
for (i in 1 : max(lst$assign)) {
  tmp2 <- colSums(lst[lst$assign==i,])<2
  tmp <- ifelse(tmp==tmp2 & tmp2==TRUE,TRUE,FALSE)
}
lst <- lst[,order(colSums(lst))]
k.pen <- lst[,tmp==TRUE]
colnames(k.pen) <- attributes(mm)$dimnames[[2]][-1]

lst <- lst[,order(colSums(lst))]

length(attributes(attributes(mf)$terms)$term.labels)



length(tmp[tmp==TRUE])


tmp <- ldply(lst, function(x) sum(x))


by(lst, lst$assign, function(x) x)

colnames(lst)!="assign"

lst[,colnames(lst)!="assign"]

?ddply

tmp <- lst[,lapply(.SD,sum),by=assign]

apply(tmp,2,)

tmp[,.SD<2]

?lapply

assign <- attributes(mm)$assign[-1][1:count]

lst$cols

#remove within factor comparisons


#grepl(":",colnames(mm)[-1])
#grep[colnames(mm)[-1]]


blah <- t(expand.grid(a=c(0,1),b1=c(0,1),b2=c(0,1),c1=c(0,1),c2=c(0,1)))[,-1]
names <- rownames(blah)

combn(names,3)

blah <- as.data.frame(t(blah))
t(blah[blah[,2]+blah[,3]<2,])

blah <- t(as.data.table(blah))


z<-outer(x,y, paste, sep="")
dim(z)<-NULL

x <- 
z
"a1" "b1" "c1" "a2" "b2" "c2" "a3" "b3" "c3"
z <-as.vector(t(z))



a <- rnorm(100)
#b <- rnorm(100)
b <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))
c <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))

y <- a + b + c + rnorm(100)

b <- factor(b)
c <- factor(c)

mod <- y~a*b*c
mm <- model.matrix(mod)

grepl(":",colnames(mm)[-1])

grep[colnames(mm)[-1]]

?grep

mf <- model.frame(mod)

attributes(attributes(mf)$terms)

n.var <- 3
n.expand <- 4
n.cont <- 1
lst <- list(NULL)
for(i in 1:3) {
  lst[[i]] <- c(0,1)  
}
lst <- t(expand.grid(lst))[,-1]
lst <- lst[,order(colSums(lst))]
rownames(lst) <- c("X1","X2.1","X2.2")
lst




model.frame(mod)
ifelse(unique(model.frame(mod)!=0)[,-1],1,0)
ifelse(unique(model.matrix(mod)!=0)[,-1],1,0)

unique(model.frame(mod)[,3:4])

attr(model.matrix(mod),"assign")

       x1 <- rnorm(100)
x2 <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))
x3 <- sample(c(0,1,2),100,replace=T,c(1/3,1/3,1/3))
y <- x1 + x2 + x1*x2 + rnorm(100)
x2 <- factor(x2)
x3 <- factor(x3)

mod <- y~x1*x2*x3

model.matrix(mod)
model.frame(mod)

x2.mm <- model.matrix(~x2)[,-1]
x2.1 <- x2.mm[,1]
x2.2 <- x2.mm[,2]

x2.1 <- sample(c(0,1),100,replace=T,c(1/2,1/2))
x2.2 <- sample(c(0,1),100,replace=T,c(1/2,1/2))
y <- x1 + 

?cat
b.c <- factor(b.c)

mod1 <- y ~ a*b.c
mod2 <- y~a*b
summary(lm(mod))
gemm(mod)
gemm(mod2)
summary(lm(y~a*b))
lm(mod1)

gemm(y~a*b.c)

mymod <- y ~ a*b
mymod.c <- y.c~a*factor(b.c)

model.matrix(mymod)
model.matrix(mymod.c)
model.matrix(~a+factor(b.c))

mf <- model.frame(mymod.c)

model.offset(mf)

str(mf)
model.extract(mf, "variables")

summary(lm(mf))

attributes(attributes(mf)$terms)$factor[-1,]

get_all_vars(mf)

summary(lm(mymod.c))

?model.frame
?get_all_vars


data(trees)
ff <- log(Volume) ~ log(Height) + log(Girth)
str(m <- model.frame(ff, trees))
mat <- model.matrix(ff, m)

dd <- data.frame(a = gl(3,4), b = gl(4,1,12))# balanced 2-way
options("contrasts")
model.matrix(~a+b,dd)
model.matrix(~a+b,dd,contrasts=list(a="contr.sum"))
model.matrix(~a+b,dd,contrasts=list(a="contr.sum",b="contr.poly"))
m.orth <- model.matrix(~a+b,dd, contrasts=list(a="contr.helmert"))
crossprod(m.orth)# m.orth is  ALMOST  orthogonal

values <- matrix(c(1,2,3,4,5,6,7,8,9), ncol=3)
betas <- matrix(c(1,2,3,1,2,3),nrow=2,byrow=T)

t(betas)


values%*%t(betas)



allComb<-function(x) {
  nx<-length(x)
  combs<-as.list(x)
  for(i in 2:length(x)) {
    indices<-combn(1:nx,i)
    for(j in 1:dim(indices)[2]) {
      xlist<-list(x[indices[,j]])
      combs<-c(combs,xlist)
    }
  }
  return(combs)
} 

allComb(c(1,2,3,4))


t(contr.treatment(3))

q <- c("","time")
r <- c("","diet1","diet2","diet3")
s <- c("","c1","c2")


mat <- t(unique(expand.grid(q,r,s)))
names <- apply(mat,2,function(x) paste(x, collapse=':'))
names <- sub('::', ':', names)  
names <- sub(':$', '', names)
names <- sub('^:', '', names)[-1]


names[order(laply(sapply(names, function(x) strsplit(x,':')), length))]



blah <- sapply(names, function(x) strsplit(x,':'))

laply(sapply(names, function(x) strsplit(x,':')), length)

blah


?sapply

grepl(':',names)




gsub(": $", "", paste(1:4, ifelse(is.na(foo), "", foo), sep = ", "))



apply(mat,2, function(x) length(strsplit(x,":")) ) 

strsplit(as.character(mat), ":"), length


1
names[,order()]

sub('[[:space:]]+$', '', str) ## white space, POSIX-style
sub('\\s+$', '', str, perl = TRUE) ## Perl-style white space



x <- c("a","b")
paste(x, collapse=":")
length(strsplit(x,':'))

mat <- ifelse(t(unique(expand.grid(a,b,c)))=="",0,1)



mat[,order(colSums(mat))]


z<-outer(a,b,c, paste, sep=":")


dim(z)<-NULL

fo
typeof(fo)  # R internal : "language"

class(fo <- y ~ x1*x2) # "formula"
attr(terms(fo),"term.labels")

environment(fo)
environment(as.formula("y ~ x"))
environment(as.formula("y ~ x", env = new.env()))


## Create a formula for a model with a large number of variables:
xnam <- paste0("x", 1:25)
(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))