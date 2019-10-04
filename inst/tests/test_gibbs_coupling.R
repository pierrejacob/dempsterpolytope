library(dempsterpolytope)
library(doParallel)
registerDoParallel(cores = detectCores()-2)
library(doRNG)
set.seed(3)
rm(list = ls())
# install.packages("lpSolveAPI")
# library(lpSolveAPI)

#
## Let's generate some data 
## number of categories
K <- 5
## number of observations
n <- 200
theta_dgp <- rexp(K)
theta_dgp <- theta_dgp / sum(theta_dgp)
## observations
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K) + rep(1, K)
print(freqX)
## (note that each count was incremented by one to ensure non zero entries, for simplicity)

rmaxcoupling <- function(k, theta_star1, theta_star2){
  x <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star1)
  pdf1_x <- 1/theta_star1[k]
  pdf2_x <- dempsterpolytope:::dunif_piktheta_cpp(x, k, theta_star2)
  if (runif(1) < (pdf2_x / pdf1_x)){
    return(list(pts = cbind(x,x), equal = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star2)
      pdf2_y <- 1/theta_star2[k]
      pdf1_y <- dempsterpolytope:::dunif_piktheta_cpp(y, k, theta_star1)
      reject <- (runif(1) < (pdf1_y/pdf2_y))
    }
    return(list(pts = cbind(x,y), equal = FALSE))
  }
}

omega <- 0.9
# rinit <- function() freqX / sum(freqX)
rinit <- function(){ x <- rexp(K); x <- x / sum(x) ; return(x)}
## rinit is distribution of theta_0
meeting_times <- function(freqX, lag, rinit, omega, max_iterations = 1e4){
  K <- length(freqX)
  same_a <- list()
  same_a_in_categoryk <- rep(FALSE, K)
  for (k in 1:K){
    if (freqX[k] > 0)
      same_a[[k]] <- rep(FALSE, freqX[k]) # indicator of each a's being identical in both chains
    else 
      same_a[[k]] <- TRUE
  }
  #########set LP 
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
  #########
  categories <- 1:K
  #### store points in barycentric coordinates
  # Achain1 <- list()
  # Achain2 <- list()
  # for (k in 1:K){
  #   if (freqX[k] > 0){
  #     Achain1[[k]] <- array(0, dim = c(freqX[k], K))
  #     Achain2[[k]] <- array(0, dim = c(freqX[k], K))
  #   } else {
  #     Achain1[[k]] <- array(0, dim = c(1, K))
  #     Achain2[[k]] <- array(0, dim = c(1, K))
  #   }
  # }
  # store constraints in barycentric coordinates
  # etas_chain <- array(0, dim = c(niterations, K, K))
  # etas_chain2 <- array(0, dim = c(niterations, K, K))
  ## initialization
  theta_01 <- rinit() 
  theta_02 <- rinit() 
  init_tmp1 <- initialize_pts(freqX, theta_01)
  pts1 <- init_tmp1$pts
  init_tmp2 <- initialize_pts(freqX, theta_02)
  pts2 <- init_tmp2$pts
  
  # store points
  # for (k in 1:K){
  #   if (freqX[k] > 0){
  #     Achain[[k]][1,,] <- pts[[k]]
  #     Achain2[[k]][1,,] <- pts2[[k]]
  #   } else {
  #     Achain[[k]][1,1,] <- rep(1/(K-1), K)
  #     Achain[[k]][1,1,k] <- 0
  #     Achain2[[k]][1,1,] <- rep(1/(K-1), K)
  #     Achain2[[k]][1,1,k] <- 0
  #   }
  # }
  etas1 <- do.call(rbind, init_tmp1$minratios)
  etas2 <- do.call(rbind, init_tmp2$minratios)
  # store constraints
  # etas_chain[1,,] <- etas1
  # etas_chain2[1,,] <- etas2
  ##### Advance first chain by 'lag' step
  iteration <- 0
  for (l in 1:lag){
    iteration <- iteration + 1
    ## do a Gibbs sweep
    # loop over categories
    for (k in categories){ if (freqX[k] > 0){
      # set Linear Program for this update
      # and find theta_star
      mat_cst_ <- mat_cst; icst <- 1
      for (j in setdiff(1:K, k)){ for (i in setdiff(1:K, j)){
        if (all(is.finite(etas1[j,]))){
          row_ <- (K+1)+icst; mat_cst_[row_,i] <- 1; mat_cst_[row_,j] <- -etas1[j,i]
        }
        icst <- icst + 1
      }}
      for (ik in 1:K) set.column(lpobject, ik, mat_cst_[,ik])
      vec_ <- rep(0, K); vec_[k] <- -1; set.objfn(lpobject, vec_)
      solve(lpobject); theta_star1 <- get.variables(lpobject)
      pts_k <- dempsterpolytope:::runif_piktheta_cpp(freqX[k], k, theta_star1)
      pts1[[k]] <- pts_k$pts
      etas1[k,] <- pts_k$minratios
    }}
  }
  
  ### Now do coupled Gibbs steps until meeting
  meeting <- Inf
  while (is.infinite(meeting) && iteration < max_iterations){
    iteration <- iteration + 1
    # loop over categories
    for (k in categories){ if (freqX[k] > 0){
      mat_cst_ <- mat_cst; icst <- 1
      for (j in setdiff(1:K, k)){ for (i in setdiff(1:K, j)){
        if (all(is.finite(etas1[j,]))){
          row_ <- (K+1)+icst; mat_cst_[row_,i] <- 1; mat_cst_[row_,j] <- -etas1[j,i]
        }
        icst <- icst + 1
      }}
      for (ik in 1:K) set.column(lpobject, ik, mat_cst_[,ik])
      vec_ <- rep(0, K); vec_[k] <- -1; set.objfn(lpobject, vec_)
      solve(lpobject); theta_star1 <- get.variables(lpobject)
      
      ## find other theta_star
      mat_cst_ <- mat_cst; icst <- 1
      for (j in setdiff(1:K, k)){ for (i in setdiff(1:K, j)){
        if (all(is.finite(etas2[j,]))){
          row_ <- (K+1)+icst; mat_cst_[row_,i] <- 1; mat_cst_[row_,j] <- -etas2[j,i]
        }
        icst <- icst + 1
      }}
      for (ik in 1:K) set.column(lpobject, ik, mat_cst_[,ik])
      vec_ <- rep(0, K); vec_[k] <- -1; set.objfn(lpobject, vec_)
      solve(lpobject); theta_star2 <- get.variables(lpobject)
      ## now we have theta_star1 and theta_star2
      ## with probability omega, do common RNG, otherwise max coupling
      u_ <- runif(1)
      if (u_ < omega){
        coupled_results_ <- dempsterpolytope:::crng_runif_piktheta_cpp(freqX[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
      } else {
        pts1_ <- matrix(NA, nrow = freqX[k], ncol = K)
        pts2_ <- matrix(NA, nrow = freqX[k], ncol = K)
        for (irow in 1:freqX[k]){
          res_ <- rmaxcoupling(k, theta_star1, theta_star2)
          pts1_[irow,] <- res_$pts[,1]
          pts2_[irow,] <- res_$pts[,2]
          same_a[[k]][irow] <- res_$equal
        }
        minratios1 <- apply(pts1_ / pts1_[,k], 2, min)
        minratios2 <- apply(pts2_ / pts2_[,k], 2, min)
        etas1[k,] <- minratios1
        etas2[k,] <- minratios2
        pts1[[k]] <- pts1_
        pts2[[k]] <- pts2_
        same_a_in_categoryk <- all(same_a[[k]])
      }
    }
    }
    if (all(same_a_in_categoryk)){
      meeting <- iteration
    }
  }
  rm(lpobject)
  return(meeting)
}


tv_upper_bound_estimates <- function(coupling_times, L, t)
{
  return(pmax(0,ceiling((coupling_times-L-t)/L)))
}

NREP <- 1e3
lag <- 100
meetings <- unlist(foreach(irep = 1:NREP) %dorng% {
  meeting_times(freqX, lag = lag, rinit, omega)
})
hist(unlist(meetings))

times_ <- 1:500
tv_upperbounds <- sapply(times_, function(t) mean(tv_upper_bound_estimates(meetings, lag, t)))

plot(times_, tv_upperbounds, type = "l")

# 
# # store points and constraints
# for (k in categories){
#   if (freqX[k] > 0){
#     Achain[[k]][iter_gibbs,,] <- pts[[k]]
#     Achain2[[k]][iter_gibbs,,] <- pts2[[k]]
#   } else {
#     Achain[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
#     Achain[[k]][iter_gibbs,1,k] <- 0
#     Achain2[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
#     Achain2[[k]][iter_gibbs,1,k] <- 0
#   }
# }
# etas_chain[iter_gibbs,,] <- etas
# etas_chain2[iter_gibbs,,] <- etas2
# }

# return(list(etas_chain = etas_chain, Achain = Achain))
# }

# distance_etas <- (sapply(1:niterations, function(index) sum(abs(etas_chain[index,,] - etas_chain2[index,,]))))
# plot(1:niterations, distance_etas, type = "l", log = "y")

# 
# param <- 1
# interval <- c(0.3, 0.4)
# intervalcvxp <- interval2polytope(K, param, interval)
# # 
# burnin <- 500
# postburn <- niterations - burnin
# contained_ <- rep(0, postburn)
# intersects_ <- rep(0, postburn)
# contained_2 <- rep(0, postburn)
# intersects_2 <- rep(0, postburn)
# for (index in ((burnin+1):niterations)){
#   cvxp <- etas2cvxpolytope(etas_chain[index,,])
#   cvxp2 <- etas2cvxpolytope(etas_chain2[index,,])
#   res_ <- compare_polytopes(cvxp, intervalcvxp)
#   res_2 <- compare_polytopes(cvxp2, intervalcvxp)
#   contained_[index-burnin] <- res_[1]
#   intersects_[index-burnin] <- res_[2]
#   contained_2[index-burnin] <- res_2[1]
#   intersects_2[index-burnin] <- res_2[2]
# }
# 
# # lower/upper probabilities, post burn-in
# cat(mean(contained_), mean(intersects_), "\n")
# cat(mean(contained_2), mean(intersects_2), "\n")
# # for the interval
# interval
# # on parameter
# param
# # equal to 
# (freqX/sum(freqX))[param]
