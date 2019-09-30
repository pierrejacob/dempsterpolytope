## This scripts reproduce the results for the linkage model

library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())
# define A and b such that theta = A phi + b where phi is in the interval (0,1)
A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)

# data 
# x = c(25, 3, 4, 7)
# x = 5 * c(125,18,20,34)

## well specified model
phi_dgp <- 0.5
theta_dgp <- A * phi_dgp + b
x <- sample(1:4, 10000, replace = TRUE, prob = theta_dgp)
x <- as.numeric(table(x))

## 1. Dirichlet DSM model.

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
## show empirical lower/upper CDF
plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red',
     main = 'Dirichlet DSM upper and lower CDF, nsim = 1e4',
     xlab = 'phi0', ylab = 'P(phi < phi0)', xlim = c(0, 1))
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)

## 2. Imprecise Dirchlet Model.
post_idm <- function(x, a, by = 1e-3){
  # trapezoid integration for idm posterior
  # x: data; a: alpha
  # function returns a vector of cdf evaluated per ``by''
  phi <- seq(0, 1, by = by)
  p_ <- (1/2+phi/4)^(x[1]+a[1]-1)*
    (1/4-phi/4)^(x[2]+x[3]+a[2]+a[3]-1)*
    (phi/4)^(x[4]+a[4]-1)
  vol_ <- (p_[-1] + head(p_, -1))*by/2
  cdf <- cumsum(vol_)/sum(vol_)
  return(c(0, cdf))
}

# a grid of alpha values to try
a_grid <- expand.grid(a1 = seq(0, 1, by = 0.02),
                      a2 = seq(0, 1, by = 0.02),
                      a4 = seq(0, 1, by = 0.02))
a_grid <- a_grid[which(rowSums(a_grid) == 1), ]
a_grid$a3 <- 0
a_grid <- as.matrix(a_grid[, c(1,2,4,3)])

# evaluate idm posterior at a grid of a values
by = 1e-4
idm_grid <- array(NA, dim = c(nrow(a_grid), 1/by+1))
for (i in 1:nrow(a_grid)){
  alpha <- a_grid[i, ]
  idm_grid[i, ] <- post_idm(x=x, a=alpha, by=by)
}

# pointwise lower and upper probabilities of IDM
lb_idm <- apply(idm_grid, 2, min)
ub_idm <- apply(idm_grid, 2, max)

plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'DDSM vs IDM upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
lines(seq(0, 1, by=by), lb_idm, col = 'blue', lty = 2)
lines(seq(0, 1, by=by), ub_idm, col = 'red', lty = 2)
legend('topleft', lty = c(1, 2), legend = c('DDSM', 'IDM'))

## 3. Simplex DSM model.

# define burnin and number of iterations
burnin <- 1e3
niterations_gibbs <- 1e4
# run Gibbs sampler for linkage model
samples_gibbs <- gibbs_sampler_linkage(niterations_gibbs, freqX = x, A = A, b = b)

# this produces random intervals for phi
phi_intervals_sdsm <- samples_gibbs$lu_chain[burnin:niterations_gibbs,]
# plot endpoints of intervals for phi
# matplot(samples_gibbs$lu_chain, type = "l")

## now get lower and upper CDF for a grid of values of phi in (0,1)
phi_grid_length <- 200
phi_grid <- seq(from = 0, to = 1, length.out = phi_grid_length)
cdf_lower_values <- rep(0, phi_grid_length)
cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], phi_intervals_sdsm)
  cdf_lower_values[iphi] <- cdf_res[1]
  cdf_upper_values[iphi] <- cdf_res[2]
}
##

# jpeg(filename = "inst/figure_linkage.jpg", width = 1000, height = 1000, quality = 100)
plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
lines(seq(0, 1, by=by), lb_idm, col = 'blue', lty = 2)
lines(seq(0, 1, by=by), ub_idm, col = 'red', lty = 2)
legend('topleft', lty = c(1, 2, 3), legend = c('DDSM', 'IDM', 'SDSM'))
lines(x = phi_grid, y = cdf_lower_values, col = 'blue', lty = 3)
lines(x = phi_grid, y = cdf_upper_values, col = 'red', lty = 3)
abline(v = 0.5)
# dev.off()

## oops
###
# head(phi_intervals_sdsm)
hist(phi_intervals_sdsm[,1], prob = TRUE, nclass = 100, xlim = c(0.4, 0.6))
hist(phi_intervals_ddsm[,1], prob = TRUE, add = TRUE, col = rgb(1,0,0,0.25), nclass = 100)
### 
hist(phi_intervals_sdsm[,2], prob = TRUE, nclass = 100, xlim = c(0.4, 0.6))
hist(phi_intervals_ddsm[,2], prob = TRUE, add = TRUE, col = rgb(1,0,0,0.25), nclass = 100)

# 
# 
# etas <- samples_gibbs$etas_chain[sample(1:niterations_gibbs, 1),,]
# K <- 4
# get_lower_upper <- function(etas, A, b){
#   upper <- 1
#   lower <- 0
#   K <- dim(etas)[1]
#   for (k1 in 1:K){
#     for (k2 in setdiff(1:K, k1)){
#       c <- etas[k1,k2]
#       denom <- A[k2] - c * A[k1]
#       if (denom >= 0){
#         upper <- min(upper, (c * b[k1] - b[k2]) / (denom))
#       } else {
#         lower <- max(lower, (c * b[k1] - b[k2]) / (denom))
#       }
#     }
#   }
#   return(c(lower, upper))
# }
# interval_ <- get_lower_upper(etas, A, b)
# midpoint <- runif(1, min = interval_[1], max = interval_[2])
# theta <- A * midpoint + b
# 
# for (k1 in 1:K){
#   for (k2 in setdiff(1:K, k1)){
#     print(theta[k1]/theta[k2] <= etas[k2,k1])
#   }
# }
# print(etas)
