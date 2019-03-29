# initialize points by defining a point in the simplex as theta_0
# then sampling a_n uniformly on pi_k(theta_0) if X_n == k
# also returns 'minratios', a list of K vectors, where the k-th vector contains min_{a in A_k} a_ell / a_k at element ell
#'@export
initialize_pts <- function(freqX, theta_0){
  K_ <- length(freqX)
  pts <- list()
  minratios <- list()
  for (k in 1:K_){
    if (freqX[k] > 0){
      tmp <- runif_piktheta_cpp(freqX[k], k, theta_0)
      pts[[k]] <- tmp$pts
      minratios[[k]] <- tmp$minratios
    } else { #?
      pts[[k]] <- NA
      minratios[[k]] <- rep(Inf, K_)
    }
  }
  return(list(pts = pts, minratios = minratios))
}
#
### redraw points in category k, using graph g
#'@export
refresh_pts_category <- function(g, k, freqX){
  K_ <- length(freqX)
  theta_star <- rep(0, K_)
  # minimum value among paths from k to ell ("eta star")
  # minimum_values <- rep(1, K_)
  # for (ell in setdiff(1:K_, k)){
  #   minimum_values[ell] <- distances(g, v = ell, to = k, mode = "out")
  # }
  notk <- setdiff(1:K_, k)
  minimum_values <- rep(1, K_)
  minimum_values[notk] <- distances(g, v = notk, to = k, mode = "out")[,1]
  theta_star <- exp(-minimum_values)
  theta_star[k] <- 1
  theta_star <- theta_star / sum(theta_star)
  pts_k <- runif_piktheta_cpp(freqX[k], k, theta_star)
  return(pts_k)
}

### Gibbs sampler
#'@export
gibbs_sampler <- function(niterations, freqX, theta_0){
  K_ <- length(freqX)
  if (missing(theta_0)){
    theta_0 <- freqX / sum(freqX)
  }
  categories <- 1:K_
  # store points in barycentric coordinates
  Achain <- list()
  for (k in 1:K_){
    if (freqX[k] > 0){
      Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K_))
    } else {
      Achain[[k]] <- array(0, dim = c(niterations, 1, K_))
    }
  }
  # store constraints in barycentric coordinates
  etas_chain <- array(0, dim = c(niterations, K_, K_))
  ## initialization
  init_tmp <- initialize_pts(freqX, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K_){
    if (freqX[k] > 0){
      Achain[[k]][1,,] <- pts[[k]]
    } else {
      Achain[[k]][1,1,] <- rep(1/(K_-1), K_)
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
        tmp <- refresh_pts_category(g, k, freqX)
        pts[[k]] <- tmp$pts
        etas[k,] <- tmp$minratios
        # refresh etas and graph
        # for (to_category in setdiff(1:K, k)){
        #   # etas[k,to_category] <- get_eta(pts, k, to_category)
        #   E(g, c(k, to_category))$weight <- log(etas[k,to_category])
        # }
        seqedges <- as.numeric(sapply(notk, function(x) c(k, x)))
        E(g, seqedges)$weight <- log(etas[k, notk])
      }
    }
    # store points and constraints
    for (k in categories){
      if (freqX[k] > 0){
        Achain[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K_-1), K_)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
  }
  # return points post-burnin
  return(list(etas_chain = etas_chain, Achain = Achain))
}
