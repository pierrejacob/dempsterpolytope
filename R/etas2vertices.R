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

#'@rdname etas2vertices
#'@title Convert matrix of eta values into polytope in H- and V-representations
#'@description
#' This function converts a KxK matrix "eta" to a polytope within simplex of dimension K.
#' The polytope describes a feasible set. 
#' On barycentric coordinate recall each
#' constraint is of the form: theta_ell / theta_k <= eta[k,l] i.e. theta_ell -
#' eta[k,l] theta_k  <= 0 and theta is in the simplex, i.e. sum_{j<K} theta_j <=
#' 1, -theta_j <= 0. The function writes the constraints as a matrix A and
#' a vector b such that A x <= b then calls a function of the package 'rcdd'
#' (mimicking what's done in the 'hitandrun' package) to obtain the coordinates of the
#' vertices of the polytope, i.e. the V-representation.
#'@param etas a matrix of values, such as produced by \code{\link{gibbs_sampler}}.
#'@return A list containing
#'\enumerate{
#'\item vertices_barcoord: the vertices of the polytope in barycentric coordinates
#'\item constr: a list corresponding to the H-representation, with "constr" for the matrix A, "rhs" for the vector b,
#' and "dir" containing the sign "<=" indicating the direction of all constraints.
#'}
#'@export
etas2vertices <- function(etas){
  K_ <- dim(etas)[1]
  categories <- 1:K_
  # the constraints are on the first K-1 coordinates
  # the first ones say that the feasible set is within the simplex
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K_-1, K_-1))
  b <- c(1, rep(0, K_-1))
  # then the extra constraints come from etas
  for (d in categories){
    for (j in setdiff(categories, d)){
      if (is.finite(etas[d,j])){ 
        # cccc (wA wB wC ... )' = 0
        ccc <- rep(0, K_)
        ccc[d] <- -etas[d, j]
        ccc[j] <- 1
        cc <- ccc - ccc[K_]
        b <- c(b, -ccc[K_])
        A <- rbind(A, matrix(cc[1:(K_-1)], nrow = 1))
      } else {
        # if eta is infinite, no constraint
      }
    }
  }
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  ## make H representation
  h <- rcdd::makeH(constr$constr, constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  return(list(vertices_barcoord = vertices_barcoord, constr = constr))
}
