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

## The following converts the constraint 'theta_k in [a,b]'
## into linear constraints in barycentric coordinate, i.e. of the form A x <= b.
## Then obtain the vertices of the corresponding polytope.
#'@export
interval2polytope <- function(K, param, interval){
  # we first encode the constraint that the values are in the simplex
  # the constraints are on the first K-1 coordinates
  A <- matrix(rep(1, K-1), ncol = K-1)
  A <- rbind(A, diag(-1, K-1, K-1))
  b <- c(1, rep(0, K-1))
  # then we add the constraint theta_k < b
  ccc <- rep(0, K)
  ccc[param] <- 1
  b <- c(b, interval[2] - ccc[K])
  cc <- ccc - ccc[K]
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  # and the constraint -theta_k < -a 
  ccc[param] <- -1
  b <- c(b, -interval[1] - ccc[K])
  cc <- ccc - ccc[K]
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  # then we get the vertices of polytope
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  # vertices_barcoord <- hitandrun::findVertices(constr)
  ## make H representation
  h <- rcdd::makeH(constr$constr, constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  return(list(vertices_barcoord = vertices_barcoord, constr = constr))
}

### the following function takes two convex polytopes,
## and evaluates whether the first is contained within the second (in which case contained = 1)
## and whether the two intersects
#'@export
compare_polytopes <- function(cvxp1, cvxp2){
  K_ <- dim(cvxp1$vertices_barcoord)[2]
  ## check whether polytope is contained, i.e. all vertices satisfy linear inequalities
  contained <- all(apply(cvxp1$vertices_barcoord[,1:(K_-1),drop=F], 1, function(v) all(cvxp2$constr$constr %*% v <= cvxp2$constr$rhs)))
  # check whether there is some intersection
  test_constr <- cvxp2$constr
  test_constr$constr <- rbind(test_constr$constr, cvxp1$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp1$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp1$constr$dir)
  ## make H representation
  h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  intersects <- (dim(v)[1] != 0)
  return(c(contained, intersects))
}

## take an array of etas, as produced by the function gibbs_sampler
## and a category (index between 1 and K)
## and compute whether the corresponding sets is contained / intersects 
## with the intervals [0,x] for x provided in 'xgrid'
#'@export
etas_to_lower_upper_cdf <- function(etas, category, xgrid){
  K_ <- dim(etas)[2]
  netas <- dim(etas)[1] 
  n_in_xgrid <- length(xgrid)
  # create polytopes corresponding to [0,x] on component "category"
  intervals_ <- list()
  for (igrid in 1:n_in_xgrid){
    x <- xgrid[igrid]
    intervals_[[igrid]] <- interval2polytope(K_, category, c(0, x))
  }
  iscontained_ <- matrix(FALSE, nrow = netas, ncol = n_in_xgrid)
  intersects_ <- matrix(FALSE, nrow = netas, ncol = n_in_xgrid)
  for (ieta in 1:netas){
    # get one particular "feasible set"
    eta <- etas[ieta,,]
    eta_cvxp <- etas2vertices(eta)
    # for elements in the grid
    comparison_ <- compare_polytopes(eta_cvxp, intervals_[[1]])
    iscontained_[ieta,1] <- comparison_[1]
    intersects_[ieta,1] <- comparison_[2]
    # savings come from the fact that if 
    # the polytope is contained in [0,x], it is also contained in [0,y] with x<y
    # and likewise for the "intersects" with relation
    if (n_in_xgrid>1){
      for (igrid in 2:n_in_xgrid){
        interval_ <- intervals_[[igrid]]
        if (iscontained_[ieta, igrid-1]){
          iscontained_[ieta, igrid] <- TRUE
        } else {
          iscontained_[ieta, igrid] <- all(apply(eta_cvxp$vertices_barcoord[,1:(K_-1),drop=F], 1, function(v) all(interval_$constr$constr %*% v <= interval_$constr$rhs)))
        }
        if (intersects_[ieta, igrid-1]){
          intersects_[ieta, igrid] <- TRUE
        } else {
          test_constr <- interval_$constr
          test_constr$constr <- rbind(test_constr$constr, eta_cvxp$constr$constr)
          test_constr$rhs <- c(test_constr$rhs, eta_cvxp$constr$rhs)
          test_constr$dir <- c(test_constr$dir, eta_cvxp$constr$dir)
          ## make H representation (H for ?)
          h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
          ## try to find V representation (V for Vendetta or Vertex?)
          v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
          intersects_[ieta, igrid] <- (dim(v)[1] != 0)
        }
      }  
    }
  }  
  return(list(iscontained = iscontained_, intersects = intersects_))
}

## same but parallelize the competition with foreach 
#'@export
etas_to_lower_upper_cdf_dopar <- function(etas, category, xgrid){
  K_ <- dim(etas)[2]
  netas <- dim(etas)[1] 
  n_in_xgrid <- length(xgrid)
  # create polytopes corresponding to [0,x] on component "category"
  intervals_ <- list()
  for (igrid in 1:n_in_xgrid){
    x <- xgrid[igrid]
    intervals_[[igrid]] <- interval2polytope(K_, category, c(0, x))
  }
  res_ <- foreach (ieta = 1:netas) %dopar% {
    eta <- etas[ieta,,]
    eta_cvxp <- etas2vertices(eta)
    # for elements in the grid
    comparison_ <- compare_polytopes(eta_cvxp, intervals_[[1]])
    iscontained_ <- rep(0, n_in_xgrid)
    intersects_ <- rep(0, n_in_xgrid)
    iscontained_[1] <- comparison_[1]
    intersects_[1] <- comparison_[2]
    # savings come from the fact that if 
    # the polytope is contained in [0,x], it is also contained in [0,y] with x<y
    # and likewise for the "intersects" with relation
    if (n_in_xgrid>1){
      for (igrid in 2:n_in_xgrid){
        interval_ <- intervals_[[igrid]]
        if (iscontained_[igrid-1]){
          iscontained_[igrid] <- TRUE
        } else {
          iscontained_[igrid] <- all(apply(eta_cvxp$vertices_barcoord[,1:(K_-1),drop=F], 1, function(v) all(interval_$constr$constr %*% v <= interval_$constr$rhs)))
        }
        if (intersects_[igrid-1]){
          intersects_[igrid] <- TRUE
        } else {
          test_constr <- interval_$constr
          test_constr$constr <- rbind(test_constr$constr, eta_cvxp$constr$constr)
          test_constr$rhs <- c(test_constr$rhs, eta_cvxp$constr$rhs)
          test_constr$dir <- c(test_constr$dir, eta_cvxp$constr$dir)
          ## make H representation (H for ?)
          h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
          ## try to find V representation (V for Vendetta or Vertex?)
          v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
          intersects_[igrid] <- (dim(v)[1] != 0)
        }
      }
    }
    c(iscontained_, intersects_)
  }  
  iscontained_ <- t(sapply(res_, function(x) x[1:n_in_xgrid]))
  intersects_ <- t(sapply(res_, function(x) x[(n_in_xgrid+1):(2*n_in_xgrid)]))
  return(list(iscontained = iscontained_, intersects = intersects_))
}
