# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

## function to perform one Gibbs update, of auxiliary variables corresponding to category k
## using the "shortest path approach" implemented in the igraph package
refresh_pts_category_graph <- function(g, k, counts){
  # number of categories
  K <- length(counts)
  # solution 
  theta_star <- rep(0, K)
  # minimum value among paths from k to ell
  notk <- setdiff(1:K, k)
  minimum_values <- rep(1, K)
  # compute shortest paths in the graph
  minimum_values[notk] <- igraph::distances(g, v = notk, to = k, mode = "out")[,1]
  # compute solution based on values of shortest paths 
  theta_star <- exp(-minimum_values)
  theta_star[k] <- 1
  theta_star <- theta_star / sum(theta_star)
  # draw points in the sub-simplex Delta_k(theta^{star,k})
  pts_k <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta_star)
  # return points
  return(pts_k)
}

### Gibbs sampler using the igraph package
gibbs_sampler_graph <- function(niterations, counts, theta_0){
  # number of categories
  K <- length(counts)
  # if theta_0 not specified, started from MLE 
  if (missing(theta_0)){
    theta_0 <- counts / sum(counts)
  }
  # categories
  categories <- 1:K
  # store points in barycentric coordinates
  Us <- list()
  for (k in 1:K){
    if (counts[k] > 0){
      Us[[k]] <- array(0, dim = c(niterations, counts[k], K))
    } else {
      Us[[k]] <- array(0, dim = c(niterations, 1, K))
    }
  }
  # store constraints in barycentric coordinates
  etas <- array(0, dim = c(niterations, K, K))
  ## initialization
  init_tmp <- initialize_pts(counts, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K){
    if (counts[k] > 0){
      Us[[k]][1,,] <- pts[[k]]
    } else {
      Us[[k]][1,1,] <- rep(1/(K-1), K)
      Us[[k]][1,1,k] <- 0
    }
  }
  etas_current <- do.call(rbind, init_tmp$minratios)
  g <- igraph::graph_from_adjacency_matrix(log(etas_current), mode = "directed", weighted = TRUE, diag = FALSE)
  # store constraints
  etas[1,,] <- etas_current
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (counts[k] > 0){
        notk <- setdiff(1:K, k)
        tmp <- refresh_pts_category_graph(g, k, counts)
        pts[[k]] <- tmp$pts
        etas_current[k,] <- tmp$minratios
        # refresh etas and graph
        seqedges <- as.numeric(sapply(notk, function(x) c(k, x)))
        E(g, seqedges)$weight <- log(etas_current[k, notk])
      }
    }
    # store points and constraints
    for (k in categories){
      if (counts[k] > 0){
        Us[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Us[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
        Us[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas[iter_gibbs,,] <- etas_current
  }
  return(list(etas = etas, Us = Us))
}

## Gibbs sampler that relies on the lpSolve library
gibbs_sampler_lp <- function(niterations, counts, theta_0){
  # number of categories
  K <- length(counts)
  # set LP 
  # precompute (K-1)*(K-1)
  Km1squared <- (K-1)*(K-1)
  # number of constraints in the LP: K+1 constraints for the simplex
  # and (K-1)*(K-1) constraints of the form theta_i / theta_j < eta_{j,i}
  nconstraints <- K + 1 + Km1squared
  # matrix encoding the constraints
  mat_cst <- matrix(0, nrow = nconstraints, ncol = K)
  mat_cst[1,] <- 1
  for (i in 1:K) mat_cst[1+i,i] <- 1
  # direction of constraints
  dir_ <- c("=", rep(">=", K), rep("<=", Km1squared))
  # right hand side of constraints
  rhs_ <- c(1, rep(0, K), rep(0, Km1squared))
  # create LP object
  lpobject <- lpSolveAPI::make.lp(nrow = nconstraints, ncol = K)
  # set right hand side and direction
  lpSolveAPI::set.rhs(lpobject, rhs_)
  lpSolveAPI::set.constr.type(lpobject, dir_)
  # now we have the basic LP set up and we will update it during the run of Gibbs  
  if (missing(theta_0)){
    theta_0 <- counts / sum(counts)
  }
  categories <- 1:K
  # store points in barycentric coordinates
  Us <- list()
  for (k in 1:K){
    if (counts[k] > 0){
      Us[[k]] <- array(0, dim = c(niterations, counts[k], K))
    } else {
      Us[[k]] <- array(0, dim = c(niterations, 1, K))
    }
  }
  # store constraints in barycentric coordinates
  etas <- array(0, dim = c(niterations, K, K))
  ## initialization
  init_tmp <- initialize_pts(counts, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K){
    if (counts[k] > 0){
      Us[[k]][1,,] <- pts[[k]]
    } else {
      Us[[k]][1,1,] <- rep(1/(K-1), K)
      Us[[k]][1,1,k] <- 0
    }
  }
  etas_current <- do.call(rbind, init_tmp$minratios)
  # store constraints
  etas[1,,] <- etas_current
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (counts[k] > 0){
        # set Linear Program for this update
        mat_cst_ <- mat_cst
        # find theta_star
        icst <- 1
        for (j in setdiff(1:K, k)){
          for (i in setdiff(1:K, j)){
            ## constraint of the form
            # theta_i - eta_{j,i} theta_j < 0 
            if (all(is.finite(etas_current[j,]))){
              row_ <- (K+1)+icst
              mat_cst_[row_,i] <- 1
              mat_cst_[row_,j] <- -etas_current[j,i]
            }
            icst <- icst + 1
          }
        }
        # set LP with current constraints
        for (ik in 1:K){
          lpSolveAPI::set.column(lpobject, ik, mat_cst_[,ik])
        }
        # solve LP
        vec_ <- rep(0, K)
        vec_[k] <- -1
        lpSolveAPI::set.objfn(lpobject, vec_)
        solve(lpobject)
        theta_star <- lpSolveAPI::get.variables(lpobject)
        # once we have theta_star, we can draw points in pi_k(theta_star)
        pts_k <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta_star)
        pts[[k]] <- pts_k$pts
        etas_current[k,] <- pts_k$minratios
      }
    }
    # store points and constraints
    for (k in categories){
      if (counts[k] > 0){
        Us[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Us[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
        Us[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas[iter_gibbs,,] <- etas_current
  }
  rm(lpobject)
  return(list(etas = etas, Us = Us))
}


gibbs_sampler_deprecated <- function(niterations, counts, theta_0){
  # return(gibbs_sampler_graph(niterations, counts, theta_0))
  return(gibbs_sampler_lp(niterations, counts, theta_0))
}


#'@rdname gibbs_sampler
#'@title Gibbs sampler for DS convex polytopes 
#'@description This is the main function of the package. 
#' It runs the proposed Gibbs sampler
#' for a desired number of iterations, for a given vector of counts.
#' It generates a convex polytope at each iteration, in the form 
#' of a matrix "eta". 
#' Below the number of categories is denoted by K,
#' and corresponds to the length of the input vector 'counts'.
#' Each category k has count N_k, possibly equal to zero.
#' The zeros are removed from the counts when performing the Gibbs iterations,
#' and "added back" using \code{\link{extend_us}}.
#'@param niterations a number of iterations to perform; 
#' each iteration is a full sweep of Gibbs updates for each category.
#' For help on choosing the number of iterations to perform,
#' see the function \code{\link{meeting_times}}.
#'@param counts a vector of non-negative integers containing 
#' the count data; its length defines K, the number of categories. It can include zeros.
#'@return A list with the following entries:
#' \itemize{
#' \item "etas": an array of dimension niterations x K x K. 
#' For iteration 'iteration', etas[iteration,,] is a K x K matrix
#' representing a convex polytope.
#' The feasible set at that iteration, is the set of all vectors
#' theta such that, theta_l/theta_k < etas[iteration,k,l]
#' for all k,l in {1,...,K}.
#' \item "Us": a list of K arrays. For iteration 'iteration', 
#' and category 'k' in {1,...,K}, if N_k = 0 then Us[[k]] is NA.
#' If N_k > 0 then Us[[k]][iteration,,]
#' is a matrix with N_k rows, and K columns, representing the auxiliary variables
#' "U"'s generated at that iteration.
#' } 
#'@examples
#' \dontrun{
#' gibbs_results <- gibbs_sampler(niterations = 5, counts = c(1,2,0,3))
#' gibbs_results$etas[5,,]
#' }
#'@export
gibbs_sampler <- function(niterations, counts){
  # number of categories
  K <- length(counts)
  # find zeros
  nonzeroK <- sum(counts > 0)
  whichnonzero <- which(counts > 0)
  whichzero <- which(counts == 0)
  # the function will perform Gibbs on non-zero categories
  # and then add the zeros back at the end
  nonzerocounts <- counts[counts>0]
  ## set LP 
  Km1squared <- (nonzeroK-1)^2
  # number of constraints in the LP: K+1 constraints for the simplex
  # and (K-1)*(K-1) constraints of the form theta_i / theta_j < eta_{j,i}
  nconstraints <- nonzeroK + 1 + Km1squared
  # get constraints that define the simplex
  lc <- simplex_linearconstraints(nonzeroK)
  mat_cst <- rbind(lc$constr, matrix(0, nrow = Km1squared, ncol = nonzeroK))
  rhs_ <- c(lc$rhs, rep(0, Km1squared))
  dir_ <- c(lc$dir, rep("<=", Km1squared))
  # create LP object
  lpobject <- lpSolveAPI::make.lp(nrow = nconstraints, ncol = nonzeroK)
  # set right hand side and direction
  lpSolveAPI::set.rhs(lpobject, rhs_)
  lpSolveAPI::set.constr.type(lpobject, dir_)
  # now we have the basic LP set up and we will update it during the run of Gibbs  
  # initialize points using random point in simplex
  theta_0 <- rexp(nonzeroK)
  theta_0 <- theta_0 / sum(theta_0)
  categories <- 1:nonzeroK
  # store auxiliary variables 
  Us <- list()
  for (k in 1:nonzeroK){
    Us[[k]] <- array(0, dim = c(niterations, nonzerocounts[k], nonzeroK))
  }
  # store eta-constraints 
  etas <- array(0, dim = c(niterations, nonzeroK, nonzeroK))
  ## initialization
  init_tmp <- initialize_pts(nonzerocounts, theta_0)
  # 'pts' is a list with 'nonzeroK' entries, 
  # each is a matrix of dim = nonzerocounts[k] x nonzeroK
  pts <- init_tmp$pts
  ## store points at initial iteration
  for (k in 1:nonzeroK){
    Us[[k]][1,,] <- pts[[k]]
  }
  # current eta-constraints
  etas_current <- do.call(rbind, init_tmp$minratios)
  # store eta-constraints
  etas[1,,] <- etas_current
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      # set Linear Program for this update
      mat_cst_ <- mat_cst
      # find theta_star
      icst <- 1
      for (j in setdiff(categories, k)){
        for (i in setdiff(categories, j)){
          ## constraint of the form
          # theta_i - eta_{j,i} theta_j < 0 
          if (all(is.finite(etas_current[j,]))){
            row_ <- (nonzeroK+1)+icst
            mat_cst_[row_,i] <- 1
            mat_cst_[row_,j] <- -etas_current[j,i]
          }
          icst <- icst + 1
        }
      }
      # set LP with current constraints
      for (ik in categories){
        lpSolveAPI::set.column(lpobject, ik, mat_cst_[,ik])
      }
      # solve LP
      objective_vec_ <- rep(0, nonzeroK)
      objective_vec_[k] <- -1
      lpSolveAPI::set.objfn(lpobject, objective_vec_)
      solve(lpobject)
      # construct theta_star
      theta_star <- lpSolveAPI::get.variables(lpobject)
      # once we have theta_star, we can draw points in pi_k(theta_star)
      pts_k <- dempsterpolytope:::runif_piktheta_cpp(nonzerocounts[k], k, theta_star)
      pts[[k]] <- pts_k$pts
      # update eta-constraints
      etas_current[k,] <- pts_k$minratios
    }
    # store auxiliary points and eta-constraints
    for (k in categories){
      Us[[k]][iter_gibbs,,] <- pts[[k]]
    }
    etas[iter_gibbs,,] <- etas_current
  }
  rm(lpobject)
  ## add empty categories if need be, and return etas and Us
  if (length(whichzero) == 0){
    return(list(Us = Us, etas = etas))
  } else {
    ## now we need to add back the empty categories
    return(extend_us(Us, whichnonzero, whichzero))
  }
}
