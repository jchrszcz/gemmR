# Test Beta Vectors
vec7 <- c(T,T,T,T,T,T,T)
vec6 <- c(T,F,T,T,T,T,T)
vec5 <- c(T,T,F,T,T,T,T)
vec4 <- c(T,T,F,F,F,F,T)
vec3 <- c(F,F,F,F,F,F,T)
vec2 <- c(F,F,F,F,T,F,F)
vec1 <- c(F,F,F,F,F,F,F)
vecs <- rbind(vec1,vec2,vec3,vec4,vec5,vec6,vec7)

# Create test factor matrix
# Would be supplied by gemm
a <- seq(1,10)
b <- seq(1,10)
c <- seq(1,10)
y <- seq(1,10)


# Run the function
kCorFact(levels,vecs)

# Full Model Version ... One step at a time
# Pass it a matrix of factor levels from model.frame, and a beta vector, and it will return te correction factor
kCorFact <- function(mf.factor, beta.vecs) {
  k <- dim(mf.factor)[2]
  factors <- log2(k+1)
  levels <- rbind(mf.factor, diag(k)[(factors+1):k,])
  return(apply(beta.vecs,1, function(x) sum(levels%*%x!=0)))
}

#levels <- t(expand.grid( rep( list( 0:1), factors))[-1,])
#levels <- levels[,order(colSums(levels))]
