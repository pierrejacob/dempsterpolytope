
### refresh etas and graph using Gibbs sampler 
#'@export
smc_refresh_category <- function(etas, g, freqX){
  K <- length(freqX)
  for (k in 1:K){
    if (freqX[k] > 0){
      notk <- setdiff(1:K, k)
      theta_star <- rep(0, K)
      # minimum value among paths from k to ell ("eta star")
      # minimum_values <- rep(1, K)
      # for (ell in setdiff(1:K, k)){
      #   minimum_values[ell] <- distances(g, v = ell, to = k, mode = "out")
      # }
      # 
      minimum_values <- rep(1, K)
      minimum_values[notk] <- distances(g, v = notk, to = k, mode = "out")[,1]
      
      theta_star <- exp(-minimum_values)
      theta_star[k] <- 1
      theta_star <- theta_star / sum(theta_star)
      pts_k <- runif_piktheta_cpp(freqX[k], k, theta_star)
      etas[k,] <- pts_k$minratios
      # # refresh etas and graph
      # for (to_category in setdiff(1:K, k)){
      #   # etas[k,to_category] <- get_eta(pts, k, to_category)
      #   E(g, c(k, to_category))$weight <- log(etas[k,to_category])
      # }
      seqedges <- as.numeric(sapply(notk, function(x) c(k, x)))
      E(g, seqedges)$weight <- log(etas[k, notk])
    }
  }
  return(list(g = g, etas = etas))
}


## SMC sampler
#'@export
SMC_sampler <- function(nparticles, X, K, essthreshold = 0.75){
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
  #
  if (nobs > 1){
    for (iobs in 2:nobs){
      cat("assimilating observation ", iobs, "\n")
      freqX_ <- tabulate(X[1:iobs], nbins = K)
      k_ <- X[iobs]
      notk_ <- setdiff(1:K, k_)
      # propagate particles
      for (iparticle in 1:nparticles){
        g <- graphs[[iparticle]]
        etas <- etas_particles[iparticle,,]
        theta_star <- rep(0, K)
        # minimum value among paths from k to ell ("eta star")
        # minimum_values <- rep(1, K)
        # for (ell in setdiff(1:K, k_)){
        #   minimum_values[ell] <- distances(g, v = ell, to = k_, mode = "out")
        # }
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
      # resampling if ESS is low
      ess[iobs] <- 1/(sum(weights^2)) / nparticles
      if (ess[iobs] < essthreshold){
        cat("resampling at step ", iobs, "\n")
        ancestors <- SSP_resampling_(nparticles, weights)
        graphs <- graphs[ancestors]
        logweights <- rep(0, nparticles)
        weights <- rep(1/nparticles, nparticles)
        etas_particles <- etas_particles[ancestors,,]
        # MCMC moves
        for (iparticle in 1:nparticles){
          res_ <- smc_refresh_category(etas_particles[iparticle,,], g = graphs[[iparticle]], freqX_)
          graphs[[iparticle]] <- res_$g
          etas_particles[iparticle,,] <- res_$etas
        }
      }
    }
  }
  return(list(etas_particles = etas_particles, normcst = normcst, weights = weights, ess = ess))
}
