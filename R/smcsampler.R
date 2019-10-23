# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
# Copyright (C) 2019 Pierre E. Jacob
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

### refresh etas and graph using Gibbs sampler 
#'@export
smc_refresh_category_graph <- function(etas, g, freqX){
  K <- length(freqX)
  for (k in 1:K){
    if (freqX[k] > 0){
      notk <- setdiff(1:K, k)
      theta_star <- rep(0, K)
      # 
      # minimum value among paths from k to ell ("eta star")
      minimum_values <- rep(1, K)
      minimum_values[notk] <- distances(g, v = notk, to = k, mode = "out")[,1]
      # theta_star is the intersection of theta_ell/theta_k = etastar[k,ell]
      theta_star <- exp(-minimum_values)
      theta_star[k] <- 1
      theta_star <- theta_star / sum(theta_star)
      pts_k <- runif_piktheta_cpp(freqX[k], k, theta_star)
      etas[k,] <- pts_k$minratios
      seqedges <- as.numeric(sapply(notk, function(x) c(k, x)))
      E(g, seqedges)$weight <- log(etas[k, notk])
    }
  }
  return(list(g = g, etas = etas))
}


## SMC sampler
#'@export
SMC_sampler_graph <- function(nparticles, X, K, essthreshold = 0.75, resamplingtimes = NULL, verbose = FALSE, h = NULL){
  nobs <- length(X)
  etas_particles <- array(dim = c(nparticles, K, K))
  graphs <- list()
  # initialization, by drawing a uniformly on simplex
  for (iparticle in 1:nparticles){
    a <- rexp(K)
    a <- a / sum(a)
    etas <- matrix(Inf, K, K)
    diag(etas) <- 1
    k_ <- X[1]
    notk_ <- setdiff(1:K, k_)
    for (j in notk_){
      etas[k_,j] <- a[j] / a[k_]
    }
    etas_particles[iparticle,,] <- etas  
    graphs[[iparticle]] <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  }
  # storage
  logweights <- rep(0, nparticles)
  incrweights <- rep(0, nparticles)
  weights <- rep(1/nparticles, nparticles)
  normcst <- rep(0, nobs)
  normcst[1] <- 0
  ess <- rep(1, nobs)
  resamplingtimes_sofar <- c()
  #
  if (nobs > 1){
    for (iobs in 2:nobs){
      if (verbose) cat("assimilating observation ", iobs, "\n")
      freqX_ <- tabulate(X[1:iobs], nbins = K)
      k_ <- X[iobs]
      notk_ <- setdiff(1:K, k_)
      # propagate particles
      for (iparticle in 1:nparticles){
        g <- graphs[[iparticle]]
        etas <- etas_particles[iparticle,,]
        theta_star <- rep(0, K)
        # minimum value among paths from k to ell ("eta star")
        minimum_values <- rep(1, K)
        minimum_values[notk_] <- distances(g, v = notk_, to = k_, mode = "out")[,1]
        
        theta_star <- exp(-minimum_values)
        theta_star[k_] <- 1
        theta_star <- theta_star / sum(theta_star)
        pts_k <- runif_piktheta_cpp(1, k_, theta_star)
        a <- pts_k$pts[1,]
        for (j in notk_){
          etas[k_,j] <- min(etas[k_,j], a[j] / a[k_])
          # E(g, c(k_, j))$weight <- log(etas[k_,j])
        }
        seqedges <- as.numeric(sapply(notk_, function(x) c(k_, x)))
        E(g, seqedges)$weight <- log(etas[k_,notk_])
        
        incrweights[iparticle] <- log(theta_star[k_])
        #
        etas_particles[iparticle,,] <- etas
        graphs[[iparticle]] <- g
      }
      # normalize weights
      maxlogw <- max(incrweights)
      incrweights <- exp(incrweights - maxlogw)
      normcst[iobs] <- maxlogw + log(sum(weights * incrweights))
      logweights <- logweights + log(incrweights)
      weights <- exp(logweights - max(logweights))
      weights <- weights/sum(weights)
      # resampling if ESS is low or if iobs is in the provided vector 'resamplingtimes'
      ess[iobs] <- 1/(sum(weights^2)) / nparticles
      if (((ess[iobs] < essthreshold) && is.null(resamplingtimes)) || (iobs %in% resamplingtimes)){
        if (verbose) cat("resampling at step ", iobs, "\n")
        resamplingtimes_sofar <- c(resamplingtimes_sofar, iobs)
        ancestors <- SSP_resampling_(nparticles, weights)
        graphs <- graphs[ancestors]
        logweights <- rep(0, nparticles)
        weights <- rep(1/nparticles, nparticles)
        etas_particles <- etas_particles[ancestors,,]
        # MCMC moves
        for (iparticle in 1:nparticles){
          res_ <- smc_refresh_category_graph(etas_particles[iparticle,,], g = graphs[[iparticle]], freqX_)
          graphs[[iparticle]] <- res_$g
          etas_particles[iparticle,,] <- res_$etas
        }
      }
    }
  }
  hestimator <- NULL
  if (!is.null(h)){
    # compute estimator based on weighted particles
    hestimator <- weights[1] * h(etas_particles[1,,])
    for (iparticle in 2:nparticles){
      hestimator <- hestimator + weights[iparticle] * h(etas_particles[iparticle,,])  
    }
  }
  return(list(etas_particles = etas_particles, normcst = normcst, weights = weights, ess = ess,
              essthreshold = essthreshold, resamplingtimes = resamplingtimes_sofar, hestimator = hestimator))
}

#'@export
SMC_sampler_lp <- function(nparticles, X, K, essthreshold = 0.75, resamplingtimes = NULL, verbose = FALSE, h = NULL){
  nobs <- length(X)
  categories <- 1:K
  etas_particles <- array(dim = c(nparticles, K, K))
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
  # initialization, by drawing a uniformly on simplex
  for (iparticle in 1:nparticles){
    a <- rexp(K)
    a <- a / sum(a)
    etas <- matrix(Inf, K, K)
    diag(etas) <- 1
    k_ <- X[1]
    notk_ <- setdiff(1:K, k_)
    for (j in notk_){
      etas[k_,j] <- a[j] / a[k_]
    }
    etas_particles[iparticle,,] <- etas  
  }
  # storage
  logweights <- rep(0, nparticles)
  incrweights <- rep(0, nparticles)
  weights <- rep(1/nparticles, nparticles)
  normcst <- rep(0, nobs)
  normcst[1] <- 0
  ess <- rep(1, nobs)
  resamplingtimes_sofar <- c()
  #
  if (nobs > 1){
    for (iobs in 2:nobs){
      if (verbose) cat("assimilating observation ", iobs, "\n")
      freqX_ <- tabulate(X[1:iobs], nbins = K) # could be updated recursively
      k_ <- X[iobs]
      notk_ <- setdiff(1:K, k_)
      # propagate particles
      for (iparticle in 1:nparticles){
        etas <- etas_particles[iparticle,,]
        # set linear program
        mat_cst_ <- mat_cst
        # find theta_star
        icst <- 1
        for (j in notk_){
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
        vec_[k_] <- -1
        set.objfn(lpobject, vec_)
        solvestatus <- solve(lpobject)
        theta_star <- get.variables(lpobject)
        pts_k <- dempsterpolytope:::runif_piktheta_cpp(1, k_, theta_star)
        a <- pts_k$pts[1,]
        for (j in notk_){
          etas[k_,j] <- min(etas[k_,j], a[j] / a[k_])
        }
        incrweights[iparticle] <- log(theta_star[k_])
        #
        etas_particles[iparticle,,] <- etas
      }
      # normalize weights
      maxlogw <- max(incrweights)
      incrweights <- exp(incrweights - maxlogw)
      normcst[iobs] <- maxlogw + log(sum(weights * incrweights))
      logweights <- logweights + log(incrweights)
      weights <- exp(logweights - max(logweights))
      weights <- weights/sum(weights)
      # resampling if ESS is low or if iobs is in the provided vector 'resamplingtimes'
      ess[iobs] <- 1/(sum(weights^2)) / nparticles
      if (((ess[iobs] < essthreshold) && is.null(resamplingtimes)) || (iobs %in% resamplingtimes)){
        if (verbose) cat("resampling at step ", iobs, "\n")
        resamplingtimes_sofar <- c(resamplingtimes_sofar, iobs)
        ancestors <- SSP_resampling_(nparticles, weights)
        logweights <- rep(0, nparticles)
        weights <- rep(1/nparticles, nparticles)
        etas_particles <- etas_particles[ancestors,,]
        # MCMC moves
        for (iparticle in 1:nparticles){
          etas <- etas_particles[iparticle,,]
          # loop over categories
          for (k in categories){
            if (freqX_[k] > 0){
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
              solve(lpobject)
              theta_star <- get.variables(lpobject)
              # once we have theta_star, we can draw points in pi_k(theta_star)
              pts_k <- dempsterpolytope:::runif_piktheta_cpp(freqX_[k], k, theta_star)
              # pts[[k]] <- pts_k$pts
              etas[k,] <- pts_k$minratios
            }
          }
          etas_particles[iparticle,,] <- etas
        }
      }
    }
  }
  hestimator <- NULL
  if (!is.null(h)){
    # compute estimator based on weighted particles
    hestimator <- weights[1] * h(etas_particles[1,,])
    for (iparticle in 2:nparticles){
      hestimator <- hestimator + weights[iparticle] * h(etas_particles[iparticle,,])  
    }
  }
  rm(lpobject)
  return(list(etas_particles = etas_particles, normcst = normcst, weights = weights, ess = ess,
              essthreshold = essthreshold, resamplingtimes = resamplingtimes_sofar, hestimator = hestimator))
}

#'@export
SMC_sampler <- function(nparticles, X, K, essthreshold = 0.75, resamplingtimes = NULL, verbose = FALSE, h = NULL){
  SMC_sampler_lp(nparticles, X, K, essthreshold, resamplingtimes, verbose, h)
}
