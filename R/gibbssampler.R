# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
# Copyright (C) 2019 Pierre E. Jacob
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

# initialize points by defining a point in the simplex as theta_0
# then sampling a_n uniformly on pi_k(theta_0) if X_n == k
# also returns 'minratios', a list of K vectors, where the k-th vector contains min_{a in A_k} a_ell / a_k at element ell

#'@rdname initialize_pts
#'@title Initialize auxiliary variables 
#'@description For a point theta in the simplex, a vector of counts, sample points 'u_n' uniformly in the subsimplex Delta_k(theta) 
#' where x_n = k.  
#'@param counts a vector of K integers representing counts 
#'@param theta a point in the simplex, represented by a K-vector of non-negative values summing to one
#'@return A list with 'pts', containing a list of K matrices of size N_k x K containing the auxiliary variables
#' and 'minratios', containing a list of K vectors of size K, where minratios[[k]] has entry ell equal to min_{n in I_k} u_{n,ell} / u_{n,k}
#' which are later used to define 'eta'. The notation follows that of the paper. 
#'@examples
#' initialize_pts(c(3,1,2), c(1/3,1/3,1/3))
#'@export
initialize_pts <- function(counts, theta){
  # number of categories
  K <- length(counts)
  # create list to store points
  pts <- list()
  # create list to store minimum ratio
  minratios <- list()
  # for each category
  for (k in 1:K){
    # if there are counts in that category
    if (counts[k] > 0){
      # generate sample in appropriate subsimplex
      tmp <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta)
      # store them and associated minimum ratios
      pts[[k]] <- tmp$pts
      minratios[[k]] <- tmp$minratios
    } else { # no count in that category
      pts[[k]] <- NA
      minratios[[k]] <- rep(Inf, K)
    }
  }
  # return points and associated minimum ratios ('eta')
  return(list(pts = pts, minratios = minratios))
}

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

#'@rdname gibbs_sampler
#'@title Gibbs sampler for Categorical inference 
#'@description This is the main function of the package. It runs the proposed Gibbs sampler
#' for a desired number of iterations, for a given vector of counts.
#' It generates a convex polytope at each iteration. Below the number of categories is denoted by K,
#' and corresponds to the length of the input vector 'counts'.
#'@param niterations a number of iterations to perform (each iteration is a full sweep of Gibbs updates)
#'@param counts a vector of non-negative integers containing the count data; its length defines K, the number of categories.
#'@param theta_0 (optional) a vector in the K-simplex used to initialize the sampler.
#' If missing, it is set to counts / sum(counts). 
#'@return A list containing 'etas', an array of dimension niterations x K x K
#' and 'Us', a list of K arrays, each of dimension niterations x N_k x K where N_k is the number of 
#' observations in category k. 
#' \enumerate{
#' \item etas: for iteration 'iteration', etas[iteration,,] is a K x K matrix.
#' The feasible set at that iteration, is the set of all vectors theta such that, theta_l/theta_k < etas[iteration,k,l]
#' for all k,l in {1,...,K}.
#' \item Us: for iteration 'iteration', and category 'k' in {1,...,K}, Us[[k]][iteration,,]
#' is made of N_k rows, where N_k is the k-th entry of counts, i.e. the number of counts in category k.
#' Each row contains one of the vectors u_{n} with n in I_k in the notation of the article.
#' } 
#'@examples
#' gibbs_results <- gibbs_sampler(niterations = 5, counts = c(1,2,3))
#' gibbs_results$etas[5,,]
#'@export
gibbs_sampler <- function(niterations, counts, theta_0){
  # return(gibbs_sampler_graph(niterations, counts, theta_0))
  return(gibbs_sampler_lp(niterations, counts, theta_0))
}

