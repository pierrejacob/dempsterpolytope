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

## function to sample u's uniformly on the simplex, 
## and then compute and return the associate etas, 
## where etas[k,ell] is the minimum of u[ell] / u[k] over the u's corresponding to category k
sample_uniform_etas <- function(X, K){
  # number of observations
  n <- length(X)
  # matrix of uniform a's in the simplex
  u <- matrix(rexp(K*n), ncol = K)
  u <- t(apply(u, 1, function(v) v / sum(v)))
  # create etas
  etas <- diag(1, K, K)
  for (k in 1:K){
    notk <- setdiff(1:K, k)
    u_k <- u[X == k,,drop=F]
    if (dim(u_k)[1] == 0){
      etas[k,notk] <- +Inf
    } else {
      for (ell in notk){
        etas[k,ell] <- min(u_k[,ell]/u_k[,k])
      }
    }
  }
  return(etas)
}

## function to check that the etas satisfy the inequalities:
## for all L, for all j_1, ..., j_L etas[j_1,j_2] * etas[j_2,j_3] * ... * etas[j_L,j_1] >= 1
## return TRUE if inequalities are satisfied, FALSE otherwise
check_cst_graph <- function(etas){
  g <- igraph::graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  negativecycle <- inherits(try(igraph::distances(g, mode = "out"), silent = TRUE), "try-error")
  return(!negativecycle)
}

#'@rdname rejection_sampler
#'@title Rejection sampler for DS convex polytopes
#'@description Implements a rejection sampler, 
#' as an alternative to the proposed Gibbs sampler. Only works
#' for very small data sets, otherwise the number of attempts required
#' to obtain a draw is prohibitively large.
#' The sampler draws N uniform variables in the simplex of dimension K,
#' and then checks whether the constraints are satisfied, specifically
#' by checking whether a certain graph contains negative cycles. 
#' The implementation uses the "igraph" package for that.
#'@param counts a vector of counts, of length K, that could contain zeros.
#'@param maxnattempts an integer indicating the number of trials to perform before giving up.
#'@return A list with the following entries:
#'\itemize{
#'\item "etas": a KxK matrix if procedure succeeded, otherwise NULL.
#'\item "nattempts": an integer indicating the number of trials.
#'} 
#'@examples 
#' \dontrun{
#' rejection_sampler(c(1,2,3))
#' rejection_sampler(c(6,0,5,0))
#' }
#'@export
rejection_sampler <- function(counts, maxnattempts = 1e5){
  X <- rep(1:length(counts), times = counts)
  K <- length(counts)
  accept <- FALSE
  etas <- NULL
  nattempts <- 0
  while ((!accept) && (nattempts < maxnattempts)){
    etas <- sample_uniform_etas(X, K)
    nattempts <- nattempts + 1
    accept <- check_cst_graph(etas)
  }
  if (nattempts == maxnattempts){
    return(list(etas = NULL, nattempts = nattempts))
  } else {
    return(list(etas = etas, nattempts = nattempts))
  }
}

