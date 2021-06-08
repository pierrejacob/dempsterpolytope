# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
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

#'@rdname sample_meeting_times
#'@title Sample meeting times associated with coupled lagged Gibbs chains
#'@description
#' Sample meeting times to monitor convergence
#' of the Gibbs sampler, as in the article
#' *Estimating Convergence of Markov chains with L-Lag Couplings*
#' by Niloy Biswas, Pierre E. Jacob, Paul Vanetti,
#' available at <https://arxiv.org/abs/1905.09971>.
#'
#' The coupled chains are generated using a mixture of coupled kernels;
#' with probability omega, a "common random numbers (CRN)" coupling is performed;
#' otherwise a (nearly) maximal coupling is used in a Gibbs sweep.
#'@param counts vector of counts; could include zeros. 
#'@param lag lag between the chains; defaults to 1.
#'@param omega probability of a coupled "CRN" step as opposed to maximal coupling step; defaults to 0.9.
#'@param max_iterations number of iterations after which to stop, in case meeting hasn't occurred; defaults to 1e5.
#'@param removezero remove zeros from counts before running coupled Gibbs sampler; defaults to TRUE
#'@return An integer representing the meeting time; or +Inf if meeting has not occurred before 'max_iterations'.
#'The meeting times can be turned into upper bounds on the total variation (TV) distance
#'between the Markov chain at some iteration and the limiting distribution,
#'using the function \code{\link{tv_upper_bound}}.
#'@examples
#'\dontrun{
#'nrep <- 100
#'meeting_times <- sapply(1:nrep, function(irep) sample_meeting_times(counts = c(3,2,0,1)))
#'tmax <- floor(max(meeting_times)*1.2)
#'ubounds <- sapply(1:tmax, function(t) tv_upper_bound(meeting_times, 1, t))
#'plot(x = 1:tmax, y = ubounds, type = 'l', xlab = "iteration", ylab = "TV upper bounds")
#'}
#'@export
sample_meeting_times <- function(counts, lag = 1, omega = 0.9, max_iterations = 1e5, removezero = TRUE){
  K <- length(counts) # number of categories
  if (removezero){
    counts <- counts[counts>0]
    K <- length(counts) # number of categories
  }
  rinit <- function(){ x = rexp(K); return(x/sum(x))}
  categories <- 1:K
  same_u_in_categoryk <- rep(FALSE, K) # indicates whether all variables in a category are identical
  same_u <- list() # indicates whether the auxiliary variables are identical across the chains
  for (k in 1:K){
    if (counts[k] > 0){
      same_u[[k]] <- rep(FALSE, counts[k]) # indicator of each a's being identical in both chains
    } else { 
      same_u[[k]] <- TRUE
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
  init_tmp1 <- initialize_pts(counts, theta_01)  
  pts1 <- init_tmp1$pts
  init_tmp2 <- initialize_pts(counts, theta_02)
  pts2 <- init_tmp2$pts
  # compute etas  
  etas1 <- do.call(rbind, init_tmp1$minratios)
  etas2 <- do.call(rbind, init_tmp2$minratios)
  ##### advance first chain by 'lag' steps
  iteration <- 0
  for (l in 1:lag){
    iteration <- iteration + 1
    ## do a Gibbs sweep, looping over the categories
    for (k in categories){ if (counts[k] > 0){
      # set Linear Program for this update and find associated theta_star
      mat_cst_ <- mat_cst; icst <- 1
      for (j in setdiff(1:K, k)){ 
        for (i in setdiff(1:K, j)){
        if (all(is.finite(etas1[j,]))){
          row_ <- (K+1)+icst; mat_cst_[row_,i] <- 1; mat_cst_[row_,j] <- -etas1[j,i]
        }
        icst <- icst + 1
      }}
      for (ik in 1:K) set.column(lpobject, ik, mat_cst_[,ik])
      vec_ <- rep(0, K); vec_[k] <- -1; set.objfn(lpobject, vec_)
      solve(lpobject); theta_star1 <- get.variables(lpobject)
      # using theta_star, re-draw auxiliary variables
      pts_k <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta_star1)
      pts1[[k]] <- pts_k$pts
      etas1[k,] <- pts_k$minratios
    }}
  }
  ### perform  coupled Gibbs steps until the two chains meet
  meeting <- Inf
  while (is.infinite(meeting) && iteration < max_iterations){
    iteration <- iteration + 1
    # loop over categories
    for (k in categories){ if (counts[k] > 0){
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
        coupled_results_ <- dempsterpolytope:::crng_runif_piktheta_cpp(counts[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
      } else {
        ## maximal coupling
        pts1_ <- matrix(NA, nrow = counts[k], ncol = K)
        pts2_ <- matrix(NA, nrow = counts[k], ncol = K)
        coupled_results_ <- dempsterpolytope:::maxcoupling_runif_piktheta_cpp(counts[k], k, theta_star1, theta_star2)
        pts1[[k]] <- coupled_results_$pts1
        etas1[k,] <- coupled_results_$minratios1
        pts2[[k]] <- coupled_results_$pts2
        etas2[k,] <- coupled_results_$minratios2
        same_u[[k]] <- coupled_results_$equal
        ## indicate whether all auxiliary variables coincide across two chains
        same_u_in_categoryk <- all(same_u[[k]])
      }
    }}
    ## if all auxiliary variables, in all categories, coincide across two chains
    if (all(same_u_in_categoryk)){
      ## then chains have met
      meeting <- iteration
    }
  }
  ## remove Linear Program object 
  rm(lpobject)
  ## return meeting
  return(meeting)
}


