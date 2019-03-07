## function to sample a's uniformly on the simplex, 
## and then compute and return the associate etas, 
## where etas[k,ell] is the minimum of a[ell] / a[k] over the a's corresponding to category k
sample_uniform_etas <- function(X, K){
  # number of observations
  n <- length(X)
  # matrix of uniform a's in the simplex
  a <- matrix(rexp(K*n), ncol = K)
  a <- t(apply(a, 1, function(v) v / sum(v)))
  # create etas
  etas <- diag(1, K, K)
  for (k in 1:K){
    notk <- setdiff(1:K, k)
    a_k <- a[X == k,,drop=F]
    for (ell in notk){
      etas[k,ell] <- min(a_k[,ell]/a_k[,k])
    }
  }
  return(etas)
}

## function to check that the etas satisfy the inequalities:
## for all L, for all j_1, ..., j_L etas[j_1,j_2] * etas[j_2,j_3] * ... * etas[j_L,j_1] >= 1
## return TRUE if inequalities are satisfied, FALSE otherwise
check_cst_graph <- function(etas){
  g <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  negativecycle <- inherits(try(distances(g, mode = "out"), silent = TRUE), "try-error")
  return(!negativecycle)
}

## rejection sampler
## sample etas uniformly and then check if constraints are satisfied
## returns accepted etas, for which constraints are satisfied, 
## and number of attempts it took 
#'@export
rejectionsampler <- function(X, K){
  accept <- FALSE
  etas <- NULL
  nattempts <- 0
  while (!accept){
    etas <- sample_uniform_etas(X, K)
    nattempts <- nattempts + 1
    accept <- check_cst_graph(etas)
  }
  return(list(etas = etas, nattempts = nattempts))
}

