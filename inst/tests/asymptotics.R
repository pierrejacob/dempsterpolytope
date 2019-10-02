## To do: obtain lower/upper CDF on the marginals

library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations
n <- 1000
# number of categories
K <- 4
categories <- 1:K
# data 
theta_dgp <- c(0.2, 0.4, 0.3, 0.1)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
freqX <- tabulate(X)
##
## encompassing triangle has three vertices
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")

###
niterations <- 5000
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
(proc.time() - pct)[3]
##

param <- 4
interval <- theta_dgp[param] + c(-0.01, +0.01)
intervalcvxp <- interval2polytope(K, param, interval)
# 
burnin <- 1000

nsubiterations <- 500
subiterations <- floor(seq(from = burnin+1, to = niterations, length.out = nsubiterations))
contained_ <- rep(0, nsubiterations)
intersects_ <- rep(0, nsubiterations)
for (index in 1:nsubiterations){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[subiterations[index],,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index] <- res_[1]
  intersects_[index] <- res_[2]
}

# lower/upper probabilities, post burn-in
cat(mean(contained_), mean(intersects_), "\n")


xgrid <- c(0, seq(from = theta_dgp[param]-0.1, to = theta_dgp[param]+0.1, length.out = 100), 1)
cdfs_ <- etas_to_lower_upper_cdf_dopar(samples_gibbs$etas_chain[subiterations,,], param, xgrid)
lowercdf <- colMeans(cdfs_$iscontained)
uppercdf <- colMeans(cdfs_$intersects)
uppercdf <- colMeans(cdfs_$iscontained)

# lower / upper 
plot(xgrid, lowercdf, type = "l")
lines(xgrid, uppercdf)
abline(v = theta_dgp[param], lty = 2)
abline(v = freqX[param]/n, lty = 3)

