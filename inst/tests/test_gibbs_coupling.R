rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
set.seed(4)

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
counts <- tabulate(X, nbins = K) + rep(1, K)
print(counts)
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
# rinit <- function() counts / sum(counts)
rinit <- function(){ x <- rexp(K); x <- x / sum(x) ; return(x)}
## rinit is distribution of theta_0
meeting_times <- function(counts, lag, rinit, omega, max_iterations = 1e4){
  K <- length(counts)
  same_a <- list()
  same_a_in_categoryk <- rep(FALSE, K)
  for (k in 1:K){
    if (counts[k] > 0)
      same_a[[k]] <- rep(FALSE, counts[k]) # indicator of each a's being identical in both chains
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
  # Us1 <- list()
  # Us2 <- list()
  # for (k in 1:K){
  #   if (counts[k] > 0){
  #     Us1[[k]] <- array(0, dim = c(counts[k], K))
  #     Us2[[k]] <- array(0, dim = c(counts[k], K))
  #   } else {
  #     Us1[[k]] <- array(0, dim = c(1, K))
  #     Us2[[k]] <- array(0, dim = c(1, K))
  #   }
  # }
  # store constraints in barycentric coordinates
  # etas <- array(0, dim = c(niterations, K, K))
  # etas2 <- array(0, dim = c(niterations, K, K))
  ## initialization
  theta_01 <- rinit() 
  theta_02 <- rinit() 
  init_tmp1 <- initialize_pts(counts, theta_01)
  pts1 <- init_tmp1$pts
  init_tmp2 <- initialize_pts(counts, theta_02)
  pts2 <- init_tmp2$pts
  
  # store points
  # for (k in 1:K){
  #   if (counts[k] > 0){
  #     Us[[k]][1,,] <- pts[[k]]
  #     Us2[[k]][1,,] <- pts2[[k]]
  #   } else {
  #     Us[[k]][1,1,] <- rep(1/(K-1), K)
  #     Us[[k]][1,1,k] <- 0
  #     Us2[[k]][1,1,] <- rep(1/(K-1), K)
  #     Us2[[k]][1,1,k] <- 0
  #   }
  # }
  etas1 <- do.call(rbind, init_tmp1$minratios)
  etas2 <- do.call(rbind, init_tmp2$minratios)
  # store constraints
  # etas[1,,] <- etas1
  # etas2[1,,] <- etas2
  ##### Advance first chain by 'lag' step
  iteration <- 0
  for (l in 1:lag){
    iteration <- iteration + 1
    ## do a Gibbs sweep
    # loop over categories
    for (k in categories){ if (counts[k] > 0){
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
      pts_k <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta_star1)
      pts1[[k]] <- pts_k$pts
      etas1[k,] <- pts_k$minratios
    }}
  }
  
  ### Now do coupled Gibbs steps until meeting
  meeting <- Inf
  while (is.infinite(meeting) && iteration < max_iterations){
    iteration <- iteration + 1
    # loop over categories
    for (k in categories){ if (counts[k] > 0){
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
        coupled_results_ <- dempsterpolytope:::crng_runif_piktheta_cpp(counts[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
      } else {
        pts1_ <- matrix(NA, nrow = counts[k], ncol = K)
        pts2_ <- matrix(NA, nrow = counts[k], ncol = K)
        for (irow in 1:counts[k]){
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
  meeting_times(counts, lag = lag, rinit, omega)
})
hist(unlist(meetings))

times_ <- 1:500
tv_upperbounds <- sapply(times_, function(t) mean(tv_upper_bound_estimates(meetings, lag, t)))

plot(times_, tv_upperbounds, type = "l")

# 
# # store points and constraints
# for (k in categories){
#   if (counts[k] > 0){
#     Us[[k]][iter_gibbs,,] <- pts[[k]]
#     Us2[[k]][iter_gibbs,,] <- pts2[[k]]
#   } else {
#     Us[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
#     Us[[k]][iter_gibbs,1,k] <- 0
#     Us2[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
#     Us2[[k]][iter_gibbs,1,k] <- 0
#   }
# }
# etas[iter_gibbs,,] <- etas
# etas2[iter_gibbs,,] <- etas2
# }

# return(list(etas = etas, Us = Us))
# }

# distance_etas <- (sapply(1:niterations, function(index) sum(abs(etas[index,,] - etas2[index,,]))))
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
#   cvxp <- etas2cvxpolytope(etas[index,,])
#   cvxp2 <- etas2cvxpolytope(etas2[index,,])
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
# (counts/sum(counts))[param]
