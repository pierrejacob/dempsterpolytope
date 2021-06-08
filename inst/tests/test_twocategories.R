## This script looks at the case of two categories 

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

# number of observations
n <- 25
# number of categories
K <- 2
categories <- 1:K
# data 
theta_dgp <- c(0.2, 0.8)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
counts <- tabulate(X)
##

# ## determine sensible burn-in using coupled chains 
# ## omega indicates the probability of doing a common random number move
# ## versus a "maximal coupling" move
# omega <- 0.9
# nrep <- 1e3
# lag <- 200
# meetingtimes <- unlist(foreach(irep = 1:nrep) %dorng% {
#   meeting_times(counts, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
# })
# 
# ## and we can obtain upper bounds on the Total Variation distance between 
# ## the chain at some step t and its stationary distribution
# niterations <- 300
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
# ## which we can plot
# plot(1:niterations, ubounds, type = "l", xlab = "iteration", ylab = "TV upper bounds")
# ## the plot confirms that it takes less than 50 steps to get close to stationarity
# ## so we could choose a burn-in value of 50 
# 
# ## Burn-in of 500 should be more than enough
burnin <- 500
niterations <- 1e5
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
(proc.time() - pct)[3]
##

dim(samples_gibbs$etas)
length(samples_gibbs$Us)
dim(samples_gibbs$Us[[1]])
dim(samples_gibbs$Us[[2]])

##
vertices_  <- etas_vertices(samples_gibbs$etas[1,,])
vertices_
##
# vertices_cart <- t(apply(vertices_, 1, function(v) barycentric2cartesian(v, v_cartesian)))
# plot(x = c(0,1), y = c(0,0), type = "l")
# abline(v = vertices_cart[,1], lty = 2)

# param <- 1
# nsubiterations <- 1e4
# subiterations <- floor(seq(from = burnin+1, to = niterations, length.out = nsubiterations))

# xgrid <- seq(from = 0, to = 1, length.out = 100)
# cdfs_ <- etas_to_lower_upper_cdf_dopar(samples_gibbs$etas[subiterations,,], param, xgrid)
# lowercdf <- colMeans(cdfs_$iscontained)
# uppercdf <- colMeans(cdfs_$intersects)

# lower / upper 
# plot(xgrid, lowercdf, type = "l")
# lines(xgrid, uppercdf)
# abline(v = theta_dgp[param], lty = 2)
# abline(v = counts[param]/n, lty = 3)

## Now manual implementation of the Gibbs sampler, as described in Algorithm 1 of the paper
gibbs_sampler_K2_with_uniforms <- function(niterations, counts){
  if (length(counts) != 2){
    stop("counts argument should be made of two integers")
  }
  nobs <- sum(counts)
  X <- c(rep(1, counts[1]), rep(2, counts[2]))
  # initialization
  intervals <- matrix(nrow = niterations, ncol = 2)
  U <- rep(0, nobs)
  # initialization
  U[X==1] <- runif(counts[1], min = 0, max = 0.5)
  U[X==2] <- runif(counts[2], min = 0.5, max = 1)
  current <- c(max(U[X==1]), min(U[X==2]))
  intervals[1,] <- current
  for (iteration in 2:niterations){
    U[X==1] <- runif(counts[1], min = 0, max = current[2])
    current[1] <- max(U[X==1])
    U[X==2] <- runif(counts[2], min = current[1], max = 1)
    current[2] <- min(U[X==2])
    intervals[iteration,] <- current  
  }
  return(intervals)
}

intervals_K2_ <- gibbs_sampler_K2_with_uniforms(niterations, counts)
intervals_K2_ <- intervals_K2_[subiterations,]
cdf_lowerupper <- function(theta, lu_chain){
  # intersect if luchain[1] < theta
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < theta))
  # contains if luchain[2] < theta
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < theta))
  return(c(cdf_lower, cdf_upper))
}

lowercdf_alt <- rep(0, length(xgrid))
uppercdf_alt <- rep(0, length(xgrid))

for (igrid in 1:length(xgrid)){
  res_ <- cdf_lowerupper(xgrid[igrid], intervals_K2_)  
  lowercdf_alt[igrid] <- res_[1]
  uppercdf_alt[igrid] <- res_[2]
}


plot(xgrid, lowercdf_alt, type = "l", ylab = "CDF")
lines(xgrid, uppercdf_alt)
# lines(xgrid, lowercdf_alt, lty = 2, col = "red")
# lines(xgrid, uppercdf_alt, lty = 2, col = "red")
# abline(v = counts[param]/n, lty = 3)


## Now another Gibbs sampler that directly samples the end point of the feasible set
## without explicitly instantiating uniforms

gibbs_sampler_K2 <- function(niterations, counts){
  if (length(counts) != 2){
    stop("counts argument should be made of two integers")
  }
  nobs <- sum(counts)
  X <- c(rep(1, counts[1]), rep(2, counts[2]))
  # initialization
  intervals <- matrix(nrow = niterations, ncol = 2)
  # initialization
  current <- c(0, 1)
  intervals[1,] <- current
  for (iteration in 2:niterations){
    current[1] <- current[2] * rbeta(1, counts[1], 1)
    current[2] <- current[1] + (1-current[1]) * rbeta(1, 1, counts[2])
    intervals[iteration,] <- current  
  }
  return(intervals)
}

intervals_K2_ <- gibbs_sampler_K2(niterations, counts)
intervals_K2_ <- intervals_K2_[subiterations,]
cdf_lowerupper <- function(theta, lu_chain){
  # intersect if luchain[1] < theta
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < theta))
  # contains if luchain[2] < theta
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < theta))
  return(c(cdf_lower, cdf_upper))
}

lowercdf_alt <- rep(0, length(xgrid))
uppercdf_alt <- rep(0, length(xgrid))

for (igrid in 1:length(xgrid)){
  res_ <- cdf_lowerupper(xgrid[igrid], intervals_K2_)  
  lowercdf_alt[igrid] <- res_[1]
  uppercdf_alt[igrid] <- res_[2]
}


plot(xgrid, lowercdf_alt, type = "l", ylab = "CDF")
lines(xgrid, uppercdf_alt)
lines(xgrid, lowercdf_alt, lty = 2, col = "red")
lines(xgrid, uppercdf_alt, lty = 2, col = "red")
abline(v = counts[param]/n, lty = 3)

## All these approaches seem to agree with each other

