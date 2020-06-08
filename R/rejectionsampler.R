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
  g <- igraph::graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  negativecycle <- inherits(try(igraph::distances(g, mode = "out"), silent = TRUE), "try-error")
  return(!negativecycle)
}

#'@rdname rejectionsampler
#'@title Rejection sampler to obtain convex polytopes for DS inference in Categorical distributions 
#'@description This is implemented to perform quick comparison on very small data sets; the method does not 
#'scale well with the number of counts. It samples auxiliary variables uniformly in the simplex,
#'and then checks whether constraints are satisfied. The implementation uses the igraph package.
#'@param counts a vector of counts, of length K
#'It returns accepted states only, and 
#' and the number of attempts it took 
#'@examples 
#' \dontrun{
#' rejectionsampler(c(1,2,3,1,2,3,2,2,1), 3)
#' }
#'@export
rejectionsampler <- function(counts){
  X <- rep(1:length(counts), times = counts)
  K <- length(counts)
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

