## This script runs coupled Gibbs sampler
## to assess convergence to stationarity
## and thus guide the choice of burn-in 

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


## number of categories
K <- 5
## number of observations
n <- 50
## data-generating parameter value
theta_dgp <- c(0.2, 0.4, 0.2, 0.1, 0.1)
## observations, simulated
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
counts <- tabulate(X, nbins = K)
print(counts)

lag <- 1
nrep <- 500
meetingtimes <- unlist(foreach(irep = 1:nrep) %dorng% {
  sample_meeting_times(counts, lag = lag)
})

summary(meetingtimes)
## most meeting occurred before 50 steps
## thus we re-try with 
lag <- 50
meetingtimes <- unlist(foreach(irep = 1:nrep) %dorng% {
  sample_meeting_times(counts, lag = lag)
})
summary(meetingtimes - lag)

## and we can obtain upper bounds on the Total Variation distance between 
## the chain at some step t and its stationary distribution
niterations <- 100
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
## which we can plot
plot(1:niterations, ubounds, type = "l", xlab = "iteration", ylab = "TV upper bounds")
## the plot confirms that it takes less than 50 steps to get close to stationarity
## so we could choose a burn-in value of 50 

