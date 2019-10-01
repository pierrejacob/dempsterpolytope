library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(2)
rm(list = ls())
K <- 4
freqX = c(25, 3, 4, 7)

##

A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)
constr <- matrix(0, nrow = 2, ncol = K)
rhs <- c(-b[1], A[1] + b[1])
constr[1,1] <- -1
constr[2,1] <- 1
dir <- c("<=", "<=")
for (j in 2:K){
  row_constr_ <- rep(0, K)
  row_constr_[1] <- -A[j]/A[1]
  row_constr_[j] <- 1
  rhs <- c(rhs, b[j] - b[1] * A[j] / A[1], -(b[j] - b[1] * A[j] / A[1]))
  constr <- rbind(constr, row_constr_, -row_constr_, deparse.level = 0)
  dir <- c(dir, "<=", "<=")
}
segment_constr <- list(constr = constr, dir = dir, rhs = rhs)

# hitandrun::findVertices(segment_constr)
# A * 0 + b
# A * 1 + b

# gibbs_sampler_lp <- function(niterations, freqX, theta_0){
niterations <- 1e4
phi_0 <- 0.5
theta_0 <- A * phi_0 + b
##

K <- length(freqX)
# set LP 
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
## add segment constraints
nconstraints <- nconstraints + nrow(segment_constr$constr)
mat_cst <- rbind(mat_cst, segment_constr$constr)
dir_ <- c(dir_, segment_constr$dir)
rhs_ <- c(rhs_, segment_constr$rhs)
# create LP object
lpobject <- make.lp(nrow = nconstraints, ncol = K)
# set right hand side and direction
set.rhs(lpobject, rhs_)
set.constr.type(lpobject, dir_)
# now we have the basic LP set up and we will update it during the run of Gibbs  
if (missing(theta_0)){
  theta_0 <- freqX / sum(freqX)
}
categories <- 1:K
# store points in barycentric coordinates
Achain <- list()
for (k in 1:K){
  if (freqX[k] > 0){
    Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K))
  } else {
    Achain[[k]] <- array(0, dim = c(niterations, 1, K))
  }
}
# store constraints in barycentric coordinates
etas_chain <- array(0, dim = c(niterations, K, K))
## initialization
init_tmp <- initialize_pts(freqX, theta_0)
pts <- init_tmp$pts
# store points
for (k in 1:K){
  if (freqX[k] > 0){
    Achain[[k]][1,,] <- pts[[k]]
  } else {
    Achain[[k]][1,1,] <- rep(1/(K-1), K)
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
      # print(lpobject)
      solve(lpobject)
      theta_star <- get.variables(lpobject)
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
      Achain[[k]][iter_gibbs,1,] <- rep(1/(K-1), K)
      Achain[[k]][iter_gibbs,1,k] <- 0
    }
  }
  etas_chain[iter_gibbs,,] <- etas
}
rm(lpobject)
# return(list(etas_chain = etas_chain, Achain = Achain))
# }

## redefine segment constraints in dimension K-1
constr <- matrix(0, nrow = 2, ncol = K-1)
rhs <- c(-b[1], A[1] + b[1])
constr[1,1] <- -1
constr[2,1] <- 1
dir <- c("<=", "<=")
for (j in 2:(K-1)){
  row_constr_ <- rep(0, K-1)
  row_constr_[1] <- -A[j]/A[1]
  row_constr_[j] <- 1
  rhs <- c(rhs, b[j] - b[1] * A[j] / A[1], -(b[j] - b[1] * A[j] / A[1]))
  constr <- rbind(constr, row_constr_, -row_constr_, deparse.level = 0)
  dir <- c(dir, "<=", "<=")
}
segment_constr <- list(constr = constr, dir = dir, rhs = rhs)

intersects_with_segment <- function(eta, segment_constr){
  cvxp <- etas2cvxpolytope(eta)
  test_constr <- segment_constr
  test_constr$constr <- rbind(test_constr$constr, cvxp$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp$constr$dir)
  ## make H representation (H for ?)
  h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
  ## try to find V representation (V for Vendetta or Vertex?)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  intersects <- (dim(v)[1] != 0)
  return(intersects)
}

burnin <- 1e3
etas <- etas_chain[burnin:niterations,,]
intersects_segment_ <- sapply(1:dim(etas)[1], function(i) intersects_with_segment(etas[i,,], segment_constr))
sum(intersects_segment_)
##
get_lu <- function(eta, segment_constr){
  cvxp <- etas2cvxpolytope(eta)
  test_constr <- segment_constr
  test_constr$constr <- rbind(test_constr$constr, cvxp$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp$constr$dir)
  ## make H representation (H for ?)
  h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
  ## try to find V representation (V for Vendetta or Vertex?)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  seg <- v[,-c(1,2)]
  # hitandrun::findVertices(test_constr)
  seg <- cbind(seg, 1-rowSums(seg))
  phis_ <- apply(seg, 1, function(v) (v[1] - b[1]) / A[1])
  return(sort(phis_))
}

nsegments <- sum(intersects_segment_)
segments_ <- matrix(nrow = nsegments, ncol = 2)
for (isegment in 1:nsegments){
  segments_[isegment,] <- get_lu(etas[isegment,,], segment_constr)
}
##


rdirichlet <- function (n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

## function to generate samples from the Dirichlet DSM (DDSM) model 
sample_ddsm <- function(x, nsim = 1e4, batch = nsim*100){
  # Rejection algorithm for Dirichlet DSM
  # x: 4 dimensional data vector
  # nsim: total samples wanted
  # batch: sample size per batch
  
  sample_phi <- function(n, data = x){
    z = rdirichlet(n, alpha = c(data, 1))
    s_ = cbind(4*pmax(z[, 1]-1/2, z[, 4]), 1 - 4*pmax(z[, 2], z[, 3]))
    return(s_[s_[,2]>s_[,1],])
  }
  
  s_ <- sample_phi(batch, data = x); ns = nrow(s_)
  if (ns == 0){
    stop('Acceptance rate too low, change batch size')
  }
  
  while(ns < nsim){
    s_ <- rbind(s_, sample_phi(batch, data = x))
    ns <- nrow(s_)
  }
  return(s_[1:nsim,])
}

## draw 1e4 feasible sets
phi_intervals_ddsm = sample_ddsm(x = freqX, nsim = 1e4)

## now get lower and upper CDF for a grid of values of phi in (0,1)
phi_grid_length <- 100
phi_grid <- seq(from = 0, to = 1, length.out = phi_grid_length)
cdf_lower_values <- rep(0, phi_grid_length)
cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], segments_)
  cdf_lower_values[iphi] <- cdf_res[1]
  cdf_upper_values[iphi] <- cdf_res[2]
}
##

# jpeg(filename = "inst/figure_linkage.jpg", width = 1000, height = 1000, quality = 100)
plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
legend('topleft', lty = c(1, 2, 3), legend = c('DDSM', 'SDSM'))
lines(x = phi_grid, y = cdf_lower_values, col = 'blue', lty = 3)
lines(x = phi_grid, y = cdf_upper_values, col = 'red', lty = 3)
# dev.off()


burnin <- 1e3
niterations_gibbs <- 1e4
# run Gibbs sampler for linkage model
samples_gibbs_linkage <- gibbs_sampler_linkage(niterations_gibbs, freqX, A = A, b = b)
# this produces random intervals for phi
phi_intervals_sdsm <- samples_gibbs_linkage$lu_chain[burnin:niterations_gibbs,]
# plot endpoints of intervals for phi
# matplot(samples_gibbs$lu_chain, type = "l")

## now get lower and upper CDF for a grid of values of phi in (0,1)
cdf_lower_values_previous <- rep(0, phi_grid_length)
cdf_upper_values_previous <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], phi_intervals_sdsm)
  cdf_lower_values_previous[iphi] <- cdf_res[1]
  cdf_upper_values_previous[iphi] <- cdf_res[2]
}
##

plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
legend('topleft', lty = c(1, 2, 3), legend = c('DDSM', 'SDSM', 'Previous'))
lines(x = phi_grid, y = cdf_lower_values, col = 'blue', lty = 2)
lines(x = phi_grid, y = cdf_upper_values, col = 'red', lty = 2)
lines(x = phi_grid, y = cdf_lower_values_previous, col = 'orange', lty = 3)
lines(x = phi_grid, y = cdf_upper_values_previous, col = 'orange', lty = 3)
