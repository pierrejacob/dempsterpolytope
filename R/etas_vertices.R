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

#'@rdname etas_vertices
#'@title Convert matrix of eta values into list of vertices
#'@description
#' This function converts a KxK matrix "eta" into
#' a matrix listing the vertices of the associated polytope within the simplex of dimension K.
#' The matrix "eta" describes
#' constraints of the form: theta_l / theta_k <= eta[k,l],
#' or equivalently theta_l -
#' eta[k,l] theta_k  <= 0, on top of the additional constraint that 
#' theta must be in the simplex: i.e. sum_{j<K} theta_j <=
#' 1, -theta_j <= 0. 
#' The function calls the package \code{rcdd}
#' (mimicking what's done in the \code{hitandrun} package) to obtain the coordinates of the
#' vertices of the polytope, i.e. the V-representation.
#'@param etas a matrix of values, such as produced by \code{\link{gibbs_sampler}}
#' or \code{\link{rejection_sampler}}.
#'@return The vertices of the polytope in barycentric coordinates, in an M x K matrix with M rows if there are M vertices.
#'@examples 
#' \dontrun{
#' etas <- rejection_sampler(c(1,0,3,2))$etas
#' etas_vertices(etas)
#' }
#'@export
etas_vertices <- function(etas){
  K_ <- dim(etas)[1]
  categories <- 1:K_
  # # the constraints are on the first K-1 coordinates
  # # the first ones say that the feasible set is within the simplex
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
  lc <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  ## make H representation
  h <- rcdd::makeH(lc$constr, lc$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  return(vertices_barcoord)
}
