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

#'@rdname initialize_pts
#'@title Initialize auxiliary variables given theta
#'@description For a point theta in the simplex, and a vector of counts (N_1,...,N_K),
#' sample for each k, N_k points uniformly in the subsimplex Delta_k(theta), where N_k > 0,
#' and compute the minimum ratio u_{n,l} / u_{n,k} for all l in {1,...,K}. 
#'@param counts a vector of K integers representing counts, possibly including zeros
#'@param theta a point in the simplex, represented by a K-vector of non-negative values summing to one
#'@return A list with the following entries.
#'\itemize{
#' \item "pts": a list of K matrices of size N_k x K containing the 
#' auxiliary variables, or NA if N_k = 0
#' \item "minratios", a list of K vectors of size K,
#' where minratios[[k]] has entry ell equal to min_{n in I_k} u_{n,ell} / u_{n,k}
#' which can be used to define 'eta_{k,ell}'. Equals +Inf if N_k = 0. 
#' }
#'@examples
#' \dontrun{
#' initialize_pts(c(3,0,2), c(1/3,1/3,1/3))
#' }
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
