library(montecarlodsm)
set_my_theme()
set.seed(2)
rm(list = ls())

## get data from the Lawrence et al paper
K <- 4
freqX <- c(25, 3, 4, 7)
# define A and b such that theta = A phi + b where phi is in the interval (0,1)
A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)

# define burnin and number of iterations
burnin <- 1e3
niterations_gibbs <- 1e4
# run Gibbs sampler for linkage model
samples_gibbs <- gibbs_sampler_linkage(niterations_gibbs, freqX, A = A, b = b)

# this produces random intervals for phi
lu_chain <- samples_gibbs$lu_chain[burnin:niterations_gibbs,]
# plot endpoints of intervals for phi
matplot(lu_chain, type = "l")

## now get lower and upper CDF for a grid of values of phi in (0,1)
phi_grid_length <- 200
phi_grid <- seq(from = 0, to = 1, length.out = phi_grid_length)
cdf_lower_values <- rep(0, phi_grid_length)
cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], lu_chain)
  cdf_lower_values[iphi] <- cdf_res[1]
  cdf_upper_values[iphi] <- cdf_res[2]
}
##
plot(phi_grid, cdf_lower_values, type = "l")
lines(phi_grid, cdf_upper_values)

## now we can compare with the Dirichlet-DSM implemented by Robin
library(gtools)
linkage_ddsm <- function(x, nsim = 1e4, batch = nsim*100){
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

phi_ddsm = linkage_ddsm(x = freqX, nsim = 1e4)
plot(ecdf(phi_ddsm[, 1]), col = 'red',
     main = 'upper and lower CDF, nsim = 1e4',
     xlab = 'phi0', ylab = 'P(phi < phi0)', xlim = c(0, 1))
plot(ecdf(phi_ddsm[, 2]), col = 'blue', add = T)

## we can now overlay both lower CDFs and upper CDFs

lines(phi_grid, cdf_lower_values, col = "blue", lty = 2)
lines(phi_grid, cdf_upper_values, col = "red", lty = 2)

## here the Simplex-DSM are in dashed lines, Dirichlet-DSM in full lines
## and blue/red is for lower/upper CDFs
