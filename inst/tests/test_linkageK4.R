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


niterations <- 5e4
warmup <- 1e3
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
##

nsubiterations <- 20000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))

## segment constraint
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
## check that it works
# hitandrun::findVertices(segment_constr)
# A * 0 + b
# A * 1 + b

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

etas <- samples_gibbs$etas_chain[subiterations,,]
## that could be parallelized...
intersects_segment_ <- sapply(1:dim(etas)[1], function(i) intersects_with_segment(etas[i,,], segment_constr))
sum(intersects_segment_)

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
  segments_[isegment,] <- get_lu(etas[which(intersects_segment_)[isegment],,], segment_constr)
}

hist(segments_[,1])
hist(segments_[,2])

### compare to other methods 
## 1. Dirichlet DSM model.
x <- freqX
# function to sample n times from Dirichlet(alpha) distribution
# (taken verbatim from the 'gtools' package)
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
phi_intervals_ddsm = sample_ddsm(x = x, nsim = 1e4)

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

## Now previous linkage SDSM
# define burnin and number of iterations
burnin <- 1e3
niterations_gibbs <- 1e4
# run Gibbs sampler for linkage model
samples_gibbs_linkage <- gibbs_sampler_linkage(niterations_gibbs, freqX = x, A = A, b = b)
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

jpeg(filename = "inst/figure_linkage.jpg", width = 1000, height = 1000, quality = 100)
plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
legend('topleft', lty = c(1, 2, 3), legend = c('DDSM', 'SDSM', 'Previous'))
lines(x = phi_grid, y = cdf_lower_values, col = 'blue', lty = 2)
lines(x = phi_grid, y = cdf_upper_values, col = 'red', lty = 2)
lines(x = phi_grid, y = cdf_lower_values_previous, col = 'orange', lty = 3)
lines(x = phi_grid, y = cdf_upper_values_previous, col = 'orange', lty = 3)
dev.off()

