p =  10    # number of variables
K = 3      # number of models/networks
n = 100    # number of samples per group

# Generate the sparsity pattern for all variables
# all precision matrices are identical in this toy example
ix = vector("list", p)
for (i in 1:p){
  ix[[i]] = matrix(1, K, p-1)
}
m <- matrix(1:p,nrow=p) %*% matrix(rep(1,p),ncol=p) - matrix(rep(1,p),nrow=p) %*% matrix(1:p,ncol=p)
S <- exp(-abs(m)) # inverse of S is sparse

# y is the vector specifying each model
x <- vector("list", K)
y <- vector("list", K)
Sigma <- vector("list",K)
for (k in 1:K){
  Sigma[[k]] <- S 
  x[[k]] <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma[[k]])
  x[[k]] <- scale(x[[k]], center = T, scale = T)
  y[[k]] <- rep(k, n)  
}
trainX <- do.call(rbind, x) # training data
trainY <- unlist(y)        

fit.joint.1 <- JSEM(trainX, trainY, ix, lambda = 5*log(p)/n)
Info <- vector("list", K)
for (k in 1:K){
  Info[[k]] = zeroInd(fit.joint.1$Adj[[k]], 1)$zeroArr
}
fit.joint.2 <- multi.glasso(trainX, trainY, lambda = rep(0.1*log(p)/n, K), zero=Info)


