## Sample meeting times to implement TV upper bounds as in https://arxiv.org/abs/1905.09971
## with probability omega, common random numbers are used in a Gibbs sweep
## otherwise a (nearly) maximal coupling is used in a Gibbs sweep
## freqX contains the data, i.e. the counts of X_1,...,X_K
## lag refers to the lag employed in the coupling
## rinit is the distribution of theta_0, used to draw the auxiliary variables at the 
## initial step
## and max_iterations is used to cut the while loop if for some reason meetings do not occur
#'@export
meeting_times <- function(freqX, lag, rinit, omega, max_iterations = 1e5){
  K <- length(freqX) # number of categories
  categories <- 1:K
  same_a_in_categoryk <- rep(FALSE, K) # indicates whether all variables in a category are identical
  same_a <- list() # indicates whether the auxiliary variables are identical across the chains
  for (k in 1:K){
    if (freqX[k] > 0){
      same_a[[k]] <- rep(FALSE, freqX[k]) # indicator of each a's being identical in both chains
    } else { 
      same_a[[k]] <- TRUE
    }
  }
  ######### setup Linear Program (LP) 
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
  # now we have the basic LP set up and we will update it during the run of the Gibbs sampler  
  ## initialization
  theta_01 <- rinit() # initial theta_0 for both chains
  theta_02 <- rinit() 
  # draw auxiliary variables in the partition defined by theta_0 within the simplex
  init_tmp1 <- initialize_pts(freqX, theta_01)  
  pts1 <- init_tmp1$pts
  init_tmp2 <- initialize_pts(freqX, theta_02)
  pts2 <- init_tmp2$pts
  # compute etas  
  etas1 <- do.call(rbind, init_tmp1$minratios)
  etas2 <- do.call(rbind, init_tmp2$minratios)
  ##### advance first chain by 'lag' steps
  iteration <- 0
  for (l in 1:lag){
    iteration <- iteration + 1
    ## do a Gibbs sweep, looping over the categories
    for (k in categories){ if (freqX[k] > 0){
      # set Linear Program for this update and find associated theta_star
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
      # using theta_star, re-draw auxiliary variables
      pts_k <- runif_piktheta_cpp(freqX[k], k, theta_star1)
      pts1[[k]] <- pts_k$pts
      etas1[k,] <- pts_k$minratios
    }}
  }
  ### perform  coupled Gibbs steps until the two chains meet
  meeting <- Inf
  while (is.infinite(meeting) && iteration < max_iterations){
    iteration <- iteration + 1
    # loop over categories
    for (k in categories){ if (freqX[k] > 0){
      ## find the two "theta_star"
      # find first theta_star
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
      # find second theta_star
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
      ## now that we have theta_star1 and theta_star2
      ## with probability omega, do Gibbs step with common RNG, 
      ## otherwise do Gibbs step with maximal coupling
      u_ <- runif(1)
      if (u_ < omega){
        ## common random numbers
        coupled_results_ <- crng_runif_piktheta_cpp(freqX[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
      } else {
        ## maximal coupling
        pts1_ <- matrix(NA, nrow = freqX[k], ncol = K)
        pts2_ <- matrix(NA, nrow = freqX[k], ncol = K)
        coupled_results_ <- maxcoupling_runif_piktheta_cpp(freqX[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
        same_a[[k]] <- coupled_results_$equal
        ## indicate whether all auxiliary variables coincide across two chains
        same_a_in_categoryk <- all(same_a[[k]])
      }
    }}
    ## if all auxiliary variables, in all categories, coincide across two chains
    if (all(same_a_in_categoryk)){
      ## then chains have met
      meeting <- iteration
    }
  }
  ## remove Linear Program object 
  rm(lpobject)
  ## return meeting
  return(meeting)
}


