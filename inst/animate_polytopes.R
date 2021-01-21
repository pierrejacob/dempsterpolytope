## this scripts produces an animation
## showing two trajectories of convex polytopes within the simplex
## contracting towards one another
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
library(gganimate)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)

K <- 3
counts <- c(2,1,3)
niterations <- 1e3
results <- gibbs_sampler(niterations, counts)
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 20))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(results$etas[iter,,])
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
g <- create_plot_triangle(graphsettings = graphsettings)
gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 0.25, alpha = .25, fill = 'black', colour = 'black')
gpolytopes

categories <- 1:K
v1 <- graphsettings$v_cartesian[[1]]
v2 <- graphsettings$v_cartesian[[2]]
v3 <- graphsettings$v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black") + theme_void()
# g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
# gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2), size = 0.25, alpha = .25, fill = 'black', colour = 'black')
# gpolytopes_anim <- gpolytopes + transition_states(iter, 4, 1) + ease_aes('sine-in-out')
# animate(gpolytopes_anim, nframes = 400, fps = 20)

## animate coupling 

coupled_gibbs_sampler <- function(niterations, counts, omega){
  K <- length(counts) # number of categories
  rinit <- function(){ x = rexp(K); return(x/sum(x))}
  categories <- 1:K
  same_a_in_categoryk <- rep(FALSE, K) # indicates whether all variables in a category are identical
  same_a <- list() # indicates whether the auxiliary variables are identical across the chains
  for (k in 1:K){
    if (counts[k] > 0){
      same_a[[k]] <- rep(FALSE, counts[k]) # indicator of each a's being identical in both chains
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
  init_tmp1 <- initialize_pts(counts, theta_01)  
  pts1 <- init_tmp1$pts
  init_tmp2 <- initialize_pts(counts, theta_02)
  pts2 <- init_tmp2$pts
  # compute etas  
  etas1 <- do.call(rbind, init_tmp1$minratios)
  etas2 <- do.call(rbind, init_tmp2$minratios)
  # store constraints
  etas1_history <- array(0, dim = c(niterations, K, K))
  etas2_history <- array(0, dim = c(niterations, K, K))
  etas1_history[1,,] <- etas1
  etas2_history[1,,] <- etas2
  ### perform  coupled Gibbs steps until the two chains meet
  for (iteration in 2:niterations){
    # loop over categories
    for (k in categories){ if (counts[k] > 0){
      ## find the two "theta_star"
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
        same_a[[k]] <- coupled_results_$equal
        ## indicate whether all auxiliary variables coincide across two chains
        same_a_in_categoryk <- all(same_a[[k]])
      }
    }}
    if (all(same_a_in_categoryk)){
      ## then chains have met
      meeting <- iteration
    }
    etas1_history[iteration,,] <- etas1
    etas2_history[iteration,,] <- etas2
  }
  ## remove Linear Program object 
  rm(lpobject)
  ## return meeting
  return(list(etas1 = etas1_history, etas2 = etas2_history))
}


niterations <- 20
cgibbsresults <- coupled_gibbs_sampler(niterations, counts, omega = 0.8)

subiter <- floor(seq(from = 1, to = niterations, length.out = 20))
cvxpolytope_cartesian.df1 <- data.frame()
cvxpolytope_cartesian.df2 <- data.frame()
for (iter in subiter){
  cvx1 <- etas2cvxpolytope(cgibbsresults$etas1[iter,,])
  cvx2 <- etas2cvxpolytope(cgibbsresults$etas2[iter,,])
  cvx_cartesian1 <- t(apply(cvx1$vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  cvx_cartesian2 <- t(apply(cvx2$vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_1 <- colMeans(cvx_cartesian1)
  o_1 <- order(apply(sweep(cvx_cartesian1, 2, average_1, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian1 <- cvx_cartesian1[o_1,]
  average_2 <- colMeans(cvx_cartesian2)
  o_2 <- order(apply(sweep(cvx_cartesian2, 2, average_2, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian2 <- cvx_cartesian2[o_2,]
  cvxpolytope_cartesian.df1 <- rbind(cvxpolytope_cartesian.df1, data.frame(cvx_cartesian1, iter = iter))
  cvxpolytope_cartesian.df2 <- rbind(cvxpolytope_cartesian.df2, data.frame(cvx_cartesian2, iter = iter))
}
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black") + theme_void()
gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df1, aes(x = X1, y = X2), size = 0.25, alpha = .25, fill = 'black', colour = 'black')
gpolytopes <- gpolytopes + geom_polygon(data = cvxpolytope_cartesian.df2, aes(x = X1, y = X2), size = 0.25, alpha = .25, fill = 'orange', colour = 'black')
gpolytopes_anim <- gpolytopes + transition_states(iter, 4, 1, wrap = F) + ease_aes('sine-in-out')
gpolytopes_anim

# animate(gpolytopes_anim, nframes = 400, fps = 20)


