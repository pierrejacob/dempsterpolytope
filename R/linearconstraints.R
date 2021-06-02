## returns the linear constraints corresponding
## to the simplex in a given dimension
#'@export
simplex_linearconstraints <- function(dimension){
  mat <- rbind(rep(1, dimension), diag(1, dimension, dimension))
  dir <- c("=", rep(">=", dimension))
  rhs <- c(1, rep(0, dimension))
  return(list(constr = mat, dir = dir, rhs = rhs))
}

## returns the linear constraints corresponding to theta_l / theta_k <= eta[k,l]
## skipping rows of etas that contain 'Inf' values
#'@export
eta_linearconstraints <- function(eta){
  ## number of rows of eta with all finite values
  n_nonvac <- sum(apply(eta, 1, function(v) all(is.finite(v))))
  nconstraints <- n_nonvac * (dim(eta)[2] - 1)
  mat <- matrix(0, nrow = nconstraints, ncol = dim(eta)[2])
  icst <- 0
  for (k in 1:(dim(eta)[1])){
    if (all(is.finite(eta[k,]))){
      for (ell in setdiff(1:dim(eta)[2], k)){
        icst <- icst + 1
        mat[icst,k] <- -eta[k,ell]
        mat[icst,ell] <- 1
      }
    }
  }
  dir <- rep("<=", nconstraints)
  rhs <- rep(0, nconstraints)
  return(list(constr = mat, dir = dir, rhs = rhs))
}

## return linear constraints obtained by combining
## two sets of linear constraints, given in the form
## of lists with keys 'constr', 'dir', 'rhs'
#'@export
concatenate_linearconstraints <- function(lc1, lc2){
  return(list(constr = rbind(lc1$constr, lc2$constr), rhs = c(lc1$rhs, lc2$rhs), dir = c(lc1$dir, lc2$dir))
  )
}