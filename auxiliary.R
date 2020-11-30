#' @param  index The node-sepcific group index
#' @details For example, 
#' index = cbind(do.call(cbind, rep(list(c(1,1,2,2)), 10)), do.call(cbind, rep(list(c(1,2,1,2)), 10)))
#' dimension of index: K by p
#' @return A list with 
#' \item{v.index}{The index indicating group membership of each variables.}
#' \item{g.index}{The index for graphs.}
#' \item{x.index}{The index for columns of the design matrix X.}
#' \item{b.index}{The index for recovering the beta coefficients to the correct order.}
##--------------------------------------------\
sort.variables <- function(
  index # the group index matrix of K by p-1
){
  K = nrow(index)
  p = ncol(index) + 1
  
  len = apply(index, 2, function(x)  length(unique(x)) ) 
  
  g.index = matrix(rep(1:K, p-1), nrow = K, ncol = p-1)
  x.index = order(c(t(do.call(rbind, rep(list(1:ncol(index)), K)))))
  
  #initialize the variable index
  v.index = index
  for (j in 2:ncol(index)) {
    v.index[, j] = v.index[, j] + cumsum(len)[j-1]
  }
  v.index = c(v.index)
  
  # re-order the variable index so that they are monotone
  new.order = order(v.index)
  v.index = v.index[new.order]
  x.index = x.index[new.order]
  g.index = g.index[new.order]
  b.index = index
  for (j in 1:ncol(index)) {
    b.index[, j] = order(order(index[, j]))
  }
  
  res = list(v.index = v.index, g.index = g.index, x.index = x.index, b.index = b.index)
  return(res)
}

#' #' Get the nearest positive definite matrix 
#' pd <- function(A, zeta=0.1){
#'   if (sum(A != t(A)) > 0){
#'     stop("This method only works for symmetric A!")
#'   }
#' 
#'   p <- dim(A)[1]
#'   diag(A) <- rep(0, p)
#'   diag(A) <- abs(min(eigen(A)$values)) + zeta
#'   Ainv <- chol2inv(chol(A))
#'   Ainv <- cov2cor(Ainv)
#'   A <- chol2inv(chol(Ainv))
#'   return(list(A = A, Ainv = Ainv))
#' }

#' Function to threshold values at the group levels.
#' @param v The coefficient to be thresholded.
#' @param g A vector indicating group membership of each variable.
#' @param thr The threshold. 
#' @details Also see the function `sort.variables`.
thr.group <- function(v, g, thr){
  # group level thresholding based on group l2 norm
  g.v = tapply(v, as.factor(g), function(x) paste(x, collapse=", ") )
  g.v = sapply(g.v, strsplit, split=", ")
  g.v = lapply(g.v, as.numeric)
  norm2 = unlist(lapply(g.v, function(x) rep(sqrt(sum(x^2)), length(x))))
  v[which(norm2 < thr)] = 0
  return(v)
}

#' Function to threshold values within groups.
#' @param v The coefficient to be thresholded.
#' @param g A vector indicating group membership of each variable.
#' @param thr The threshold. 
thr.coord <- function(v, g, thr){
  # coordinate level thresholding when there is group misspecification
  g.v = tapply(v, as.factor(g), function(x) paste(x, collapse=", ") )
  g.v = sapply(g.v, strsplit, split=", ")
  g.v = lapply(g.v, as.numeric)
  direction = lapply(g.v, function(x) {ifelse(x==0, 0, x/sqrt(sum(x^2)))})  
  direction = unlist(direction)
  names(direction) = NULL
  v[which(abs(direction)<thr)] = 0
  return(v)
}

#' Symmetrize a matrix
#' @param A A matrix
#' @details We take the average of `A` and its transpose to symmetrize a matrix. 
symmetrize <- function(A, eps = 1e-06){
  A <- (A + t(A))/2
  A[abs(A)<eps] <- 0
  diag(A) <- 0
  return(A)
}

#' Trace of a matrix
matTr <- function(z) sum(diag(z))


#' Obtain the indices of entries of the inverse covariance to be constrained to be zero.
#' @param Amat The adjacency matrix of p by p
#' @param r A scalar between 0 and 1, indicating whether to obtain all the indices 
#' for which entries are zero.
#' @param eps Tolerance for calling an edge in the inverse covariance. Default is 1e-06.
#' @return A list with entries
#' \item{zeroArr}{The indices for which we will zero out. Format of this variable is the same as 
#' that used in `which(a, arr.ind = TRUE)`.}
#' \item{zeroMat}{The known 0's as a matrix.}
#' \item{oneMat}{The known 1's as a matrix.}
zeroInd <- function(Amat, r, eps=1e-06){
  if (!isSymmetric(Amat)){
    stop("This method only works for symmetric matrix!")
  }
  p <- dim(Amat)[1]
  oneMat <- matrix(0, p, p)
  zeroMat <- matrix(0, p, p)
  
  one.pos <- which(abs(Amat)>=eps, arr.ind = TRUE)
  zero.pos <- which(abs(Amat)<eps, arr.ind = TRUE)
  
  zero.pos <- zero.pos[which(zero.pos[,1] > zero.pos[,2]) ,]
  sel.zero <- sample(seq(1, dim(zero.pos)[1]), r * dim(zero.pos)[1], replace = FALSE) 
  zeroMat[zero.pos[sel.zero, ]] <- 1
  zeroMat <- zeroMat + t(zeroMat)  
  zeroArr <- zero.pos[sel.zero, ]
  
  out <- list()
  out$zeroArr = zeroArr
  out$zeroMat = zeroMat
  
  if (dim(one.pos)[1] == 0){
    warning("The matrix is zero!")
    out$oneMat = matrix(0, p, p)
  } else 
  {
    one.pos <- one.pos[which(one.pos[,1] > one.pos[,2]) ,]
    if (is.null(dim(one.pos))){
      one.pos = matrix(one.pos, nrow = 1)
    }
    
    sel.one <- sample(seq(1, dim(one.pos)[1]), r * dim(one.pos)[1], replace = FALSE) 
    oneMat[one.pos[sel.one, ]] <- 1
    oneMat <- oneMat + t(oneMat)
    diag(oneMat) <- 0
    
    out$oneMat = oneMat  
  }
  
  return(out)  
}
#' Function to estimate multiple adjacency matrices using graphical lasso.
#' @param trainX An n by p matrix of training data with rows ordered by group labels
#' @param trainY A length n vector indicating group labels of input data `trainX`
#' @param lambda A scalar as the penalty parameter in each glasso problem
#' @param zero (Optional) indices of entries of inverse covariance to be constrained to be zero.
#' Since we estimate multiple graphs, this should be a list of length n, one for each covariance.
#' Default is NULL.
#' @param BIC Whether to calculate the bic.score.
#' @param eps Tolerance for calling an edge in the inverse covariance. Default is 1e-06.
#' @details 
#' This function is equivalent to separate graphical lasso for each group when `zero=NULL`. 
#' It can be used after applying `JSEM` to refit the precision matrices based on the structural
#' constraints. 
#' @return A list with entries
#' \item{Omega}{A list of estimated precision matrices.}
#' \item{Adj}{A list of adjacency matrices corresponding to `Omega`.}
#' \item{BIC}{The BIC score corresponding to the input `lambda`.}
#' \item{lambda}{The input lambda parameter.}
multi.glasso <- function(
  trainX,    
  trainY,      
  lambda,     
  zero = NULL, 
  BIC = FALSE,
  eps = 1e-06
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  n = as.numeric(table(trainY))
  
  #penalty needed for glasso
  if (length(lambda)==K) {rho = lambda} else { 
    rho = rep(lambda, K)}
  
  #Initialize the estimated precision and adjacency matrix
  Omega.hat = vector("list", K)
  Ahat = vector("list", K)
  
  #Whether there are entries that need to be constrained to zero
  if (is.null(zero)){
    zero = rep(list(zero), K)
  }
  
  if (max(sapply(zero, length)) == p*(p-1)){
    stop("One or more matrices are constrained to be zero")
  }  
  
  bic.score = rep(0, K)
  
  for (k in 1:K) {
    Ahat[[k]] = matrix(0, p, p)
    data <- trainX[which(trainY == k), ]
    
    empcov <- cov(data) #empirical cov
    while (kappa(empcov) > 1e+2){
      empcov = empcov + 0.05 * diag(p)
    }
    
    fit <- glasso(empcov, rho = rho[k], zero = zero[[k]], penalize.diagonal=FALSE, maxit = 30)
    
    Omega.hat[[k]] = (fit$wi + t(fit$wi))/2
    Ahat[[k]][abs(Omega.hat[[k]])>eps] = 1
    diag(Ahat[[k]]) = 0
    
    if (BIC){
      bic.score[k] = matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]])) + log(n[k]) * sum(Ahat[[k]])/(2*n[k])
    }
    
  }
  
  out = list(Omega = Omega.hat, 
             Adj = Ahat, 
             BIC = bic.score, 
             lambda = lambda)
  return(out)
}

#' Select the tuning parameter for JSEM
#' 
#' @param trainX An n by p matrix of training data with rows ordered by group labels.
#' @param testX An nt by p matrix of test data with rows ordered by group labels. 
#' @param trainY A length n vector indicating group labels of input data `trainX`.
#' @param testY A length nt vector indicating group labels of input data `testX`.
#' @param index  A K by p matrix of node-specific grouping structure
#' @param lambda A scalar as the penalty parameter in each glasso problem
#' @return A list with
#' \item{BIC}{The BIC scores for the input lambda sequence.}
#' \item{likelihood}{The log-likelihood for the input lambda sequence.}
sel.lambda.jsem <- function(
  trainX,
  testX, 
  trainY, 
  testY, 
  index, 
  lambda 
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  Info <- vector("list", K)
  
  # Find the sample size for each category
  n <- rep(0, K)
  for (k in 1:K){
    n[k] <- nrow(trainX[which(trainY == k), ]) 
  }
  
  N <- length(lambda)
  bic.score <- rep(0, N)
  likelihood <- rep(0, N)
  for (j in 1:N){
    cat("The ", j, "-th step in tuning... \n")
    Ahat <- JSEM(trainX, trainY, index, lambda=lambda[j])$Adj
    for (k in 1:K){
      Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
    }
    
    #In the second step of Joint group lasso, we use 0.1*log(p)/n as the default penalty for each graph.
    fit <- multi.glasso(trainX, trainY, lambda = 0.1*log(p)/n, zero = Info, BIC = T)
    bic.score[j] <- sum(fit$BIC)
    
    for (k in 1:K){
      data <- testX[which(testY == k), ]      
      empcov <- cov(data) 
      while (kappa(empcov) > 1e+2){
        empcov = empcov + 0.05 * diag(p)      
      }   
      likelihood[j] = likelihood[j] + matTr(empcov %*% fit$Omega[[k]]) - log(det(fit$Omega[[k]]))
    }
  }
  
  out <- list(BIC = bic.score, likelihood = likelihood)
  return(out)
}


#'Function to perform selection of the estimated network based on stability criteria
#' @references Meinshausen and Buhlmann - Stability selection - JRSSB - 2010
#' @param X An n by p matrix of training data with rows ordered by group labels.
#' @param Y A length n vector indicating group labels of input data `X`.
#' @param cnt The number of subsampling.
#' @param lastar The oracle lambda found in graphical lasso
#' @return 
#' \item{mat}{A list of selected matrices}
#' \item{count}{The actual number of replicates tried, which might be smaller than `cnt` if 
#' error occurs in some runs.}
stabsel.jsem <- function(X, Y, index, cnt, lastar) {
  p = ncol(X)
  K <- length(unique(Y))
  
  # translate the data matrix into a list
  trainX <- vector("list",K)
  for (k in 1:K){
    trainX[[k]] <- X[Y==k,]
  }
  X <- trainX
  n = lapply(X, nrow)
  
  X1 = vector("list",K)
  X2 = vector("list", K)
  sel.mat = vector("list", K)
  for (k in 1:K){
    sel.mat[[k]] = matrix(0, p, p)
  }
  count = 0
  for (i in 1:cnt) {
    model.1 = NULL 
    model.2 = NULL 
    for (k in 1:K){
      ind.1 = sample(seq(1, n[[k]]), n[[k]]/2, F)
      ind.2 = seq(1, n[[k]])[match(seq(1, n[[k]]), ind.1, 0) == 0]
      X1[[k]] = X[[k]][ind.1, ]
      X2[[k]] = X[[k]][ind.2, ]
      model.1 = c(model.1, rep(k, length(ind.1)))
      model.2 = c(model.2, rep(k, length(ind.2)))
    }
    tmp.1 = try(JSEM(trainX=do.call(rbind, X1), trainY=model.1, index, lambda=lastar))
    tmp.2 = try(JSEM(trainX=do.call(rbind, X2), trainY=model.2, index, lambda=lastar))
    
    if (inherits(tmp.1, "try-error") || inherits(tmp.2, "try-error")){
      warning("There might be some error!")
      next;
    }
    
    for (k in 1:K){
      sel.mat[[k]] = sel.mat[[k]] + tmp.1$Adj[[k]] + tmp.2$Adj[[k]]
    }
    
    count = count + 1
  }
  
  return(list(mat = sel.mat, count = count))
}

#' Joint estimation of multiple graphical models
#' @references Guo, Jian, et al. "Joint estimation of multiple graphical models." Biometrika 98.1 (2011): 1-15.
#' @param  trainX data
#' @param  trainY labels for categories (1, 2, 3,...)
#' @param  lambda_value: tuning parameter
#' @details This is the code from Guo et al. (2011). Make sure loading "glasso" package before using the code.
#' @return A list with 
#' \item{Omega}{The estimated precision matrices.}
#' \item{S}{The empirical covariance matrices.}
#' \item{lambda}{The input tuning parameter.}
CGM_AHP_train <- function(
  trainX,
  trainY,
  lambda_value,
  adaptive_weight = array(1, c(length(unique(trainY)), ncol(trainX), ncol(trainX))),
  limkappa = 1e+6  #limit for condition number of the sample cov
){
  ## Set the general paramters
  K <- length(unique(trainY))
  p <- ncol(trainX)
  diff_value <- 1e+10
  count <- 0
  tol_value <- 0.01
  max_iter <- 30
  
  ## Set the optimizaiton parameters
  OMEGA <- array(0, c(K, p, p))
  S <- array(0, c(K, p, p))
  OMEGA_new <- array(0, c(K, p, p))
  nk <- rep(0, K)
  
  ## Initialize Omega
  for (k in seq(1, K)) {
    idx <- which(trainY == k)
    S[k, , ] <- cov(trainX[idx, ])
    if (kappa(S[k, , ]) > limkappa) {
      S[k, , ] <- S[k, , ] + 0.01 * diag(p)
    }
    tmp <- solve(S[k, , ])
    OMEGA[k, , ] <- tmp
    nk[k] <- length(idx)
  }
  
  ## Start loop
  while ((count < max_iter) & (diff_value > tol_value)) {
    tmp <- apply(abs(OMEGA), c(2, 3), sum)
    tmp[abs(tmp) < 1e-10] <- 1e-10
    V <- 1/sqrt(tmp)
    
    for (k in seq(1, K)) {
      penalty_matrix <- lambda_value * adaptive_weight[k, , ] * V
      obj_glasso <- glasso(S[k, , ], penalty_matrix, maxit = 100)
      OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi))/2
      #OMEGA_new[k, , ] <- obj_glasso$wi
    }
    
    ## Check the convergence
    diff_value <- sum(abs(OMEGA_new - OMEGA))/sum(abs(OMEGA))
    count <- count + 1
    OMEGA <- OMEGA_new
    #cat(count, ', diff_value=', diff_value, '\n')
  }
  
  ## Filter the noise
  list.Omega <- NULL
  for (k in seq(1, K)) {
    ome <- OMEGA[k, , ]
    ww <- diag(ome)
    ww[abs(ww) < 1e-10] <- 1e-10
    ww <- diag(1/sqrt(ww))
    tmp <- ww %*% ome %*% ww
    ome[abs(tmp) < 1e-06] <- 0
    OMEGA[k, , ] <- ome
    list.Omega[[k]] <- ome
  }
  
  output <- list()
  output$OMEGA <- list.Omega
  output$S <- S
  output$lambda <- lambda_value
  
  return(output)
}

