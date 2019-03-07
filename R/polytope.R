# convert "eta" constraints to polytope within simplex of dimension K
# the polytope describes the feasible set, on barycentric coordinate
# recall each constraint is of the form: theta_ell / theta_k <= eta[k,l]
# i.e. theta_ell - eta[k,l] theta_k  <= 0
# and theta is in the simplex, i.e. sum_{j<K} theta_j <= 1, -theta_j <= 0
## the following code writes the constraints as a matrix A and a vector b
## such that A x <= b
## then calls a function of the package 'hitandrun' to obtain
## the coordinates of the vertices of the polytope
#'@export
etas2cvxpolytope <- function(etas){
  K_ <- dim(etas)[1]
  categories <- 1:K_
  # the constraints are on the first K-1 coordinates
  # the first ones say that the feasible set is within the simplex
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K-1, K-1))
  b <- c(1, rep(0, K-1))
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
  # get vertices of polytope
  vertices_barcoord <- hitandrun::findVertices(constr)
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  return(list(vertices_barcoord = vertices_barcoord, constr = constr))
}


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
  vertices_barcoord <- hitandrun::findVertices(constr)
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  return(list(vertices_barcoord = vertices_barcoord, constr = constr))
}

### the following function takes two convex polytopes,
## and evaluates whether the first is contained within the second (in which case contained = 1)
## and whether the two intersects
#'@export
compare_polytopes <- function(cvxp1, cvxp2){
  ## check whether polytope is contained, i.e. all vertices satisfy linear inequalities
  contained <- all(apply(cvxp1$vertices_barcoord[,1:(K-1)], 1, function(v) all(cvxp2$constr$constr %*% v <= cvxp2$constr$rhs)))
  # check whether there is some intersection
  test_constr <- cvxp2$constr
  test_constr$constr <- rbind(test_constr$constr, cvxp1$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp1$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp1$constr$dir)
  test_vertices <- try(hitandrun::findVertices(test_constr), silent = T)
  intersects <- FALSE
  if (inherits(test_vertices, "try-error")){
    intersects <- FALSE
  } else {
    intersects <- TRUE
  }
  return(c(contained, intersects))
}

