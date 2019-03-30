library(montecarlodsm)
set.seed(1)
rm(list = ls())
# install.packages("lpSolveAPI")
library(lpSolveAPI)

#
## Let's generate some data 
## number of categories
K <- 10
## number of observations
n <- 50
theta_dgp <- rexp(K)
theta_dgp <- theta_dgp / sum(theta_dgp)
## observations
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K) + rep(1, K)
print(freqX)
## (note that each count was incremented by one to ensure non zero entries, for simplicity)

## Gibbs sampler that relies on the lpSolve library
## instead of igraph
gibbs_sampler_with_lpsolve <- function(niterations, freqX, theta_0){
  K_ <- length(freqX)
  # set LP 
  # precompute (K-1)*(K-1)
  Km1squared <- (K_-1)*(K_-1)
  # number of constraints in the LP: K+1 constraints for the simplex
  # and (K-1)*(K-1) constraints of the form theta_i / theta_j < eta_{j,i}
  nconstraints <- K_ + 1 + Km1squared
  # matrix encoding the constraints
  mat_cst <- matrix(0, nrow = nconstraints, ncol = K_)
  mat_cst[1,] <- 1
  for (i in 1:K_) mat_cst[1+i,i] <- 1
  # direction of constraints
  dir_ <- c("=", rep(">=", K_), rep("<=", Km1squared))
  # right hand side of constraints
  rhs_ <- c(1, rep(0, K_), rep(0, Km1squared))
  # create LP object
  directlp <- make.lp(nrow = nconstraints, ncol = K)
  # set right hand side and direction
  set.rhs(directlp, rhs_)
  set.constr.type(directlp, dir_)
  # now we have the basic LP set up and we will update it during the run of Gibbs  
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
            row_ <- (K+1)+icst
            mat_cst_[row_,i] <- 1
            mat_cst_[row_,j] <- -etas[j,i]
            icst <- icst + 1
          }
        }
        # set LP with current constraints
        for (ik in 1:K_){
          set.column(directlp, ik, mat_cst_[,ik])
        }
        # solve LP
        vec_ <- rep(0, K_)
        vec_[k] <- -1
        set.objfn(directlp, vec_)
        # print(directlp)
        solve(directlp)
        theta_star <- get.variables(directlp)
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
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K_-1), K_)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
  }
  rm(directlp)
  return(list(etas_chain = etas_chain, Achain = Achain))
}



niterations_gibbs <- 200

## run Gibbs sampler with graph tools
set.seed(1)
samples_gibbs <- gibbs_sampler(niterations_gibbs, freqX)
etas <- samples_gibbs$etas_chain[niterations_gibbs,,]
etas
## run Gibbs sampler with lpSolve
set.seed(1)
samples_gibbs_lpsolve <- gibbs_sampler_with_lpsolve(niterations_gibbs, freqX)
etas_lpsolve <- samples_gibbs_lpsolve$etas_chain[niterations_gibbs,,]
etas_lpsolve
#
summary(as.numeric(abs(etas - etas_lpsolve)))

## 
# library(microbenchmark)
# microbenchmark(samples_gibbs <- gibbs_sampler(niterations_gibbs, freqX),
#                samples_gibbs_lpsolve <- gibbs_sampler_with_lpsolve(niterations_gibbs, freqX),
#                times = 10)
# 


