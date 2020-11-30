library(glasso)
library(grpreg)

#' Joint structural estimation of multiple graphs
#' @details 
#' This function uses joint group lasso to estimate multiple graphical models, 
#' given the known grouping structure. It uses thresholding to correct for group misspecification 
#' @param trainX An n by p data matrix with rows ordered by group labels
#' @param trainY A length n vector indicating group labels of input data `trainX`
#' @param index  A K by p matrix of node-specific grouping structure
#' @param lambda The penalty parameter. Can be a scalar or a length p vector 
#' @param delta1 The thresholding parameter at the group level. Default is NULL.
#' @param delta2 The thresholding parameter within each group. Default is NULL.
#' @param eps Tolerance for calling an edge in the inverse covariance. Default is 1e-06.
#' @return A list with entries 
#' \item{Adj}{The list of estimated adjacency matrices.}
#' \item{lambda}{The lambda sequence used in estimation.}
#' 
#' }
JSEM <- function(
  trainX,       
  trainY,        
  index,       
  lambda,    
  delta1 = NULL,
  delta2 = NULL,
  eps = 1e-06
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  
  lambda.grpreg = rep(lambda, p)
  
  # list observed data by model
  design = vector("list", K)
  Ahat = vector("list", K)
  for (k in 1:K){
    design[[k]] = trainX[which(trainY == k), ]
    Ahat[[k]] = matrix(0, p, p)
  }
  
  ## Start the loop for each node
  for (i in 1:p) {
    #cat("i= ", i, "\n")
    ## Create a block diagonal matrix for the variables
    list.x = vector("list", K)
    list.y = vector("list", K)
    for (k in 1:K) {
      list.x[[k]] = design[[k]][, -i]
      list.y[[k]] = design[[k]][, i]
    }
    #Duplicate columns for each k and add these extra columns to the end of the rearranged design X;
    #In the meantime, add the variable indices to the list of all variables.
    #Need to reorganize the estimate in the end.
    X = as.matrix(bdiag(list.x))
    new_X = X
    Y = unlist(list.y)
    
    ## To get the index label 
    myindex = sort.variables(index[[i]])
    X = X[, myindex$x.index]
    
    fit = grpreg(X, Y, myindex$v.index, family = "gaussian", penalty = "grLasso", lambda = lambda.grpreg[i])
    coeff = fit$beta[-1, ]    
    
    if (!is.null(delta1)){
      coeff = thr.group(coeff, myindex$v.index, delta1)
    }
    if (!is.null(delta2)){
      coeff = thr.coord(coeff, myindex$v.index, delta2)
    }
    
    tmp = matrix(coeff, nrow = K)
    for (j in 1:dim(tmp)[2]){
      tmp[, j] = tmp[myindex$b.index[, j], j]
    }
    
    for (k in 1:K){
      Ahat[[k]][i, -i] = tmp[k,]
    }     
  }
  
  # Symmetrize the thresholded coefficients; 
  # Get the symmetric adjacency matrices; 
  Theta = lapply(Ahat, FUN = symmetrize)
  Ahat = lapply(Theta, function(u) {u[abs(u)>eps] = 1; return(u)})
  
  return(list(Adj = Ahat, 
              lambda = lambda.grpreg))
}  


