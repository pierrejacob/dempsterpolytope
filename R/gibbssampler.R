# initialize points by defining a point in the simplex as theta_0
# then sampling a_n uniformly on pi_k(theta_0) if X_n == k
# also returns 'minratios', a list of K vectors, where the k-th vector contains min_{a in A_k} a_ell / a_k at element ell
#'@export
initialize_pts <- function(freqX, theta_0){
  K <- length(freqX)
  pts <- list()
  minratios <- list()
  for (k in 1:K){
    if (freqX[k] > 0){
      tmp <- runif_piktheta_cpp(freqX[k], k, theta_0)
      pts[[k]] <- tmp$pts
      minratios[[k]] <- tmp$minratios
    } else { #?
      pts[[k]] <- NA
      minratios[[k]] <- rep(Inf, K)
    }
  }
  return(list(pts = pts, minratios = minratios))
}
#
### redraw points in category k, using graph g
#'@export
refresh_pts_category_graph <- function(g, k, freqX){
  K <- length(freqX)
  theta_star <- rep(0, K)
  # minimum value among paths from k to ell ("eta star")
  notk <- setdiff(1:K, k)
  minimum_values <- rep(1, K)
  minimum_values[notk] <- distances(g, v = notk, to = k, mode = "out")[,1]
  theta_star <- exp(-minimum_values)
  theta_star[k] <- 1
  theta_star <- theta_star / sum(theta_star)
  pts_k <- runif_piktheta_cpp(freqX[k], k, theta_star)
  return(pts_k)
}

### Gibbs sampler
#'@export
gibbs_sampler_graph <- function(niterations, freqX, theta_0){
  K <- length(freqX)
  if (missing(theta_0)){
    theta_0 <- freqX / sum(freqX)
  }
  categories <- 1:K
  # store points in barycentric coordinates
  Achain <- list()
  for (k in 1:K){
    if (freqX[k] > 0){
      Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K))
    } else {
      Achain[[k]] <- array(0, dim = c(niterations, 1, K))
    }
  }
  # store constraints in barycentric coordinates
  etas_chain <- array(0, dim = c(niterations, K, K))
  ## initialization
  init_tmp <- initialize_pts(freqX, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K){
    if (freqX[k] > 0){
      Achain[[k]][1,,] <- pts[[k]]
    } else {
      Achain[[k]][1,1,] <- rep(1/(K-1), K)
      Achain[[k]][1,1,k] <- 0
    }
  }
  etas <- do.call(rbind, init_tmp$minratios)
  g <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  # store constraints
  etas_chain[1,,] <- etas
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (freqX[k] > 0){
        notk <- setdiff(1:K, k)
        tmp <- refresh_pts_category_graph(g, k, freqX)
        pts[[k]] <- tmp$pts
        etas[k,] <- tmp$minratios
        # refresh etas and graph
        seqedges <- as.numeric(sapply(notk, function(x) c(k, x)))
        E(g, seqedges)$weight <- log(etas[k, notk])
      }
    }
    # store points and constraints
    for (k in categories){
      if (freqX[k] > 0){
        Achain[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
  }
  return(list(etas_chain = etas_chain, Achain = Achain))
}

## Gibbs sampler that relies on the lpSolve library
## instead of igraph
#'@export
gibbs_sampler_lp <- function(niterations, freqX, theta_0){
  K <- length(freqX)
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
  lpobject <- make.lp(nrow = nconstraints, ncol = K)
  # set right hand side and direction
  set.rhs(lpobject, rhs_)
  set.constr.type(lpobject, dir_)
  # now we have the basic LP set up and we will update it during the run of Gibbs  
  if (missing(theta_0)){
    theta_0 <- freqX / sum(freqX)
  }
  categories <- 1:K
  # store points in barycentric coordinates
  Achain <- list()
  for (k in 1:K){
    if (freqX[k] > 0){
      Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K))
    } else {
      Achain[[k]] <- array(0, dim = c(niterations, 1, K))
    }
  }
  # store constraints in barycentric coordinates
  etas_chain <- array(0, dim = c(niterations, K, K))
  ## initialization
  init_tmp <- initialize_pts(freqX, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K){
    if (freqX[k] > 0){
      Achain[[k]][1,,] <- pts[[k]]
    } else {
      Achain[[k]][1,1,] <- rep(1/(K-1), K)
      Achain[[k]][1,1,k] <- 0
    }
  }
  etas <- do.call(rbind, init_tmp$minratios)
  # store constraints
  etas_chain[1,,] <- etas
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (freqX[k] > 0){
        # set Linear Program for this update
        mat_cst_ <- mat_cst
        # find theta_star
        icst <- 1
        for (j in setdiff(1:K, k)){
          for (i in setdiff(1:K, j)){
            ## constraint of the form
            # theta_i - eta_{j,i} theta_j < 0 
            if (all(is.finite(etas[j,]))){
              row_ <- (K+1)+icst
              mat_cst_[row_,i] <- 1
              mat_cst_[row_,j] <- -etas[j,i]
            }
            icst <- icst + 1
          }
        }
        # set LP with current constraints
        for (ik in 1:K){
          set.column(lpobject, ik, mat_cst_[,ik])
        }
        # solve LP
        vec_ <- rep(0, K)
        vec_[k] <- -1
        set.objfn(lpobject, vec_)
        # print(lpobject)
        solve(lpobject)
        theta_star <- get.variables(lpobject)
        # once we have theta_star, we can draw points in pi_k(theta_star)
        pts_k <- montecarlodsm:::runif_piktheta_cpp(freqX[k], k, theta_star)
        pts[[k]] <- pts_k$pts
        etas[k,] <- pts_k$minratios
      }
    }
    # store points and constraints
    for (k in categories){
      if (freqX[k] > 0){
        Achain[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
  }
  rm(lpobject)
  return(list(etas_chain = etas_chain, Achain = Achain))
}

#'@export
gibbs_sampler <- function(niterations, freqX, theta_0){
  # return(gibbs_sampler_graph(niterations, freqX, theta_0))
  return(gibbs_sampler_lp(niterations, freqX, theta_0))
}

