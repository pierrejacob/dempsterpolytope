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
dempsterpolytope::rejectionsampler(counts)

## repeatedly
nsamples_rs <- 5e2
samples_rs <- foreach(irep = 1:nsamples_rs) %dorng% {
  dempsterpolytope::rejectionsampler(counts)
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


# ## check normalizing constant against SMC
# nparticles <- 2^12
# smc_res <- SMC_sampler(nparticles, X, K)
# ## log normalizing constant
# exp(sum(smc_res$normcst))
# ## seems to agree

##
niterations_gibbs <- 5100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations_gibbs, counts = counts)
elapsed <- (proc.time() - pct)[3]
elapsed
burnin <- 100

# 
param <- 2
nsubiterations_gibbs <- 500
subiterations <- floor(seq(from = burnin, to = niterations_gibbs, length.out = nsubiterations_gibbs))

xgrid <- seq(from = 0, to = 1, length.out = 100)
etas <- samples_gibbs$etas[subiterations,,]
cdf_ <- etas_to_lower_upper_cdf_dopar(etas, param, xgrid)
lowercdf <- colMeans(cdf_$iscontained)
uppercdf <- colMeans(cdf_$intersects)
plot(xgrid, lowercdf, type = "l")
lines(xgrid, uppercdf)

library(plyr)
res <- laply(samples_rs, function(l) as.matrix(l$etas))

cdf_rs <- etas_to_lower_upper_cdf_dopar(res, param, xgrid)
lowercdf_rs <- colMeans(cdf_rs$iscontained)
uppercdf_rs <- colMeans(cdf_rs$intersects)
#
plot(xgrid, lowercdf, type = "l")
lines(xgrid, uppercdf)
lines(xgrid, lowercdf_rs, lty = 2)
lines(xgrid, uppercdf_rs, lty = 2)


dim(cdf_rs$iscontained)

## and with SMC... 
# nparticles <- 2^10
# smc_res <- SMC_sampler(nparticles, X, K, essthreshold = 0.9)
# etas_particles <- smc_res$etas_particles
# weights <- smc_res$weights
# cdf_smc <- etas_to_lower_upper_cdf_dopar(etas_particles, param, xgrid)
# 
# dim(cdf_smc$iscontained)
# lowercdf_smc <- apply(cdf_smc$iscontained, 2, function(v) sum(v * weights))
# uppercdf_smc <- apply(cdf_smc$intersects, 2, function(v) sum(v * weights))
# lines(xgrid, lowercdf_smc, lty = 3)
# lines(xgrid, uppercdf_smc, lty = 3)

