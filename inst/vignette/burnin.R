## This script runs coupled Gibbs sampler
## to assess convergence to stationarity
## and thus guide the choice of burn-in 

library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

## number of categories
K <- 5
## number of observations
n <- 50
## data-generating parameter value
theta_dgp <- c(0.2, 0.4, 0.2, 0.1, 0.1)
## observations, simulated
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K)
print(freqX)

## omega indicates the probability of doing a common random number move
## versus a "maximal coupling" move
omega <- 0.9
lag <- 1
nrep <- 500
meetingtimes <- unlist(foreach(irep = 1:nrep) %dorng% {
  meeting_times(freqX, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
})

summary(meetingtimes)
## most meeting occurred before 50 steps
## thus we re-try with 
lag <- 50
meetingtimes <- unlist(foreach(irep = 1:nrep) %dorng% {
  meeting_times(freqX, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
})
summary(meetingtimes)

## and we can obtain upper bounds on the Total Variation distance between 
## the chain at some step t and its stationary distribution
niterations <- 100
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
## which we can plot
plot(1:niterations, ubounds, type = "l", xlab = "iteration", ylab = "TV upper bounds")
## the plot confirms that it takes less than 50 steps to get close to stationarity
## so we could choose a burn-in value of 50 

