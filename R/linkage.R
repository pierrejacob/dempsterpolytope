## convenience functions
## Given K-vectors A and b 
## and a KxK matrix etas
## find phi_l, phi_u
## such that for all phi in [phi_l,phi_u]
## and theta = A phi + b
## we have theta_i/theta_j <= eta[j,i] for all i \neq j 

get_lower_upper <- function(etas, A, b){
  upper <- 1
  lower <- 0
  K <- dim(etas)[1]
  for (k1 in 1:K){
    for (k2 in setdiff(1:K, k1)){
      c <- etas[k1,k2]
      denom <- A[k2] - c * A[k1]
      if (denom >= 0){
        upper <- min(upper, (c * b[k1] - b[k2]) / (denom))
      } else {
        lower <- max(lower, (c * b[k1] - b[k2]) / (denom))
      }
    }
  }
  return(c(lower, upper))
}

## same thing, but ignoring etas[k,.]
get_lower_upper_updatek <- function(etas, A, b, k){
  upper <- 1
  lower <- 0
  K <- dim(etas)[1]
  for (k1 in setdiff(1:K, k)){
    for (k2 in setdiff(1:K, k1)){
      c <- etas[k1,k2]
      denom <- A[k2] - c * A[k1]
      if (denom >= 0){
        upper <- min(upper, (c * b[k1] - b[k2]) / (denom))
      } else {
        lower <- max(lower, (c * b[k1] - b[k2]) / (denom))
      }
    }
  }
  return(c(lower, upper))
}

#'@export
gibbs_sampler_linkage <- function(niterations, freqX, phi_0, A, b){
  K_ <- length(freqX)
  if (missing(phi_0)){
    phi_0 <- 0.5
  }
  theta_0 <- A * phi_0 + b
  categories <- 1:K_
  # store points in barycentric coordinates
  Achain <- list()
  for (k in categories){
    if (freqX[k] > 0){
      Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K_))
    } else {
      Achain[[k]] <- array(0, dim = c(niterations, 1, K_))
    }
  }
  # store constraints in barycentric coordinates
  etas_chain <- array(0, dim = c(niterations, K_, K_))
  lu_chain <- matrix(0, nrow = niterations, ncol = 2)
  ## initialization
  init_tmp <- initialize_pts(freqX, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in categories){
    if (freqX[k] > 0){
      Achain[[k]][1,,] <- pts[[k]]
    } else {
      Achain[[k]][1,1,] <- rep(1/(K_-1), K_)
      Achain[[k]][1,1,k] <- 0
    }
  }
  etas <- do.call(rbind, init_tmp$minratios)
  # store constraints
  etas_chain[1,,] <- etas
  lu_chain[1,] <- get_lower_upper(etas, A, b)
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (freqX[k] > 0){
        lu_k <- get_lower_upper_updatek(etas, A, b, k)
        # if ((A * lu_k[1] + b)[k] < (A * lu_k[2] + b)[k]){
        if (A[k] > 0){
          # then max theta_k correspond to phi = upper bound
          theta_star <- A * lu_k[2] + b
        } else {
          # then max theta_k correspond to phi = lower bound
          theta_star <- A * lu_k[1] + b
        }
        ##
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
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K_-1), K_)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
    lu_chain[iter_gibbs,] <- get_lower_upper(etas, A, b)
  }
  # return points post-burnin
  return(list(etas_chain = etas_chain, Achain = Achain, lu_chain = lu_chain))
}

# for a given phi_0, and assertion {phi < phi_0}
# upper probability is proportion of time random interval [phi_lower, phi_upper] intersects with [0, phi_0]
# lower probability is proportion of time random interval [phi_lower, phi_upper] is contained in [0, phi_0]
#'@export
linkage_cdf_lowerupper <- function(phi_0, lu_chain){
  # intersect if phi_lower < phi_0
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < phi_0))
  # contains if phi_upper < phi_0 
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < phi_0))
  return(c(cdf_lower, cdf_upper))
}

