## This script looks at whether Dempster's approach
## appears to be consistent with large data sets 

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
n <- 500
# number of categories
K <- 4
categories <- 1:K
# data 
theta_dgp <- c(0.2, 0.4, 0.3, 0.1)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
counts <- tabulate(X)
##

###
niterations <- 5000
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
(proc.time() - pct)[3]
##
burnin <- 1000
nsubiterations <- 500
subiterations <- floor(seq(from = burnin+1, to = niterations, length.out = nsubiterations))
etas <- samples_gibbs$etas[subiterations,,]

param <- 4
interval <- theta_dgp[param] + c(-0.01, +0.01)

objvec <- rep(0, K)
objvec[param] <- 1

min_param <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  lpsolve_over_eta(etas[ieta,,], objvec)
}
max_param <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  -lpsolve_over_eta(etas[ieta,,], -objvec)
}

plot(min_param, max_param)
abline(v = interval)
abline(h = interval)

## intersection 
intersect_ <- !((max_param < interval[1]) | (min_param > interval[2]))
## inclusion 
contained_ <- ((min_param > interval[1]) & (max_param < interval[2]))

# lower/upper probabilities, post burn-in
cat(mean(contained_), mean(intersect_), "\n")

xgrid <- c(0, seq(from = theta_dgp[param]-0.1, to = theta_dgp[param]+0.1, length.out = 100), 1)

# lower / upper 
plot(xgrid, (ecdf(max_param))(xgrid), type = "l")
lines(xgrid, (ecdf(min_param))(xgrid), col = 'red')
abline(v = theta_dgp[param], lty = 2)
abline(v = counts[param]/n, lty = 3)

