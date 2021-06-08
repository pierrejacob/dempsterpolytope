## Rejection sampling algorithm 
## and comparison with Gibbs sampler
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
# (don't set it to be too large, because rejection sampler would become very slow; values of <= 5,6 are OK)
n <- 6
# number of categories
K <- 3
categories <- 1:K
# data 
# counts <- c(70,150,0)
# theta_dgp <- c(0.3, 0.3, 0.4)
# X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
X <- c(categories, sample(x = categories, size = n - K, replace = TRUE))
counts <- tabulate(X)
counts
##

## Rejection sampling
dempsterpolytope::rejection_sampler(counts)

## repeatedly
nsamples_rs <- 5e3
samples_rs <- foreach(irep = 1:nsamples_rs) %dorng% {
  dempsterpolytope::rejection_sampler(counts, maxnattempts = 1e5)
}

## to estimate the volume of accept region
## we can use the fact that the number of attempts before a success
## is geometric 
hist(sapply(samples_rs, function(x) x$nattempts-1))
## so an unbiased estimator of 1/Z is 
mean(sapply(samples_rs, function(x) x$nattempts))
## and an unbiased estimator of Z is 
1/(nsamples_rs/(nsamples_rs-1) * mean(sapply(samples_rs, function(x) x$nattempts-1)) + 1)
## the MLE is 
1/(mean(sapply(samples_rs, function(x) x$nattempts-1)) + 1)

##
niterations_gibbs <- 20100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations_gibbs, counts = counts)
elapsed <- (proc.time() - pct)[3]
elapsed
burnin <- 100

# 
param <- 2
nsubiterations_gibbs <- 5000
subiterations <- floor(seq(from = burnin, to = niterations_gibbs, length.out = nsubiterations_gibbs))

xgrid <- seq(from = 0, to = 1, length.out = 100)
etas <- samples_gibbs$etas[subiterations,,]

objvec <- rep(0, K)
objvec[param] <- 1

etas_min <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  lpsolve_over_eta(etas[ieta,,], objvec)
}
etas_max <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  -lpsolve_over_eta(etas[ieta,,], -objvec)
}

ecdf_lower <- ecdf(etas_max)
ecdf_upper <- ecdf(etas_min)

plot(xgrid, ecdf_lower(xgrid), type = "l")
lines(xgrid, ecdf_upper(xgrid), lty = 2)

## check against rejection sampler

etas_min_rs <- foreach(ieta = 1:nsamples_rs, .combine = c) %dopar% {
  lpsolve_over_eta(samples_rs[[ieta]]$etas, objvec)
}
etas_max_rs <- foreach(ieta = 1:nsamples_rs, .combine = c) %dopar% {
  -lpsolve_over_eta(samples_rs[[ieta]]$etas, -objvec)
}

ecdf_lower_rs <- ecdf(etas_max_rs)
ecdf_upper_rs <- ecdf(etas_min_rs)

lines(xgrid, ecdf_lower_rs(xgrid), lty = 1, col = 'red')
lines(xgrid, ecdf_upper_rs(xgrid), lty = 2, col = 'red')


