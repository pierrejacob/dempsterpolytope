## In this script, we compute lower and upper CDF of marginals 
## from the output of a Gibbs sampler for categorical distribution with arbitrary K
###
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())
## encompassing triangle has three vertices
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")

# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
theta_dgp <- c(0.5, 0.2, 0.3)
# theta_dgp <- c(0.1, 0.2, 0.3, 0.4)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
counts <- tabulate(X, nbins = K)
cat("Data:", counts, "\n")
#
niterations <- 3000
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
names(samples_gibbs)
dim(samples_gibbs$etas)
# length(samples_gibbs$Us)
# dim(samples_gibbs$Us[[1]])
#
warmup <- 2000
nsubiterations <- 1000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))

#
param <- 1

ngrid01 <- 75
grid01 <- seq(from = 0, to = 1, length.out = ngrid01)
# grid01 <- c(0, seq(from = 0.4, to = 0.6, length.out = ngrid01-2), 1)

etas <- samples_gibbs$etas[subiterations,,]

res_ <- etas_to_lower_upper_cdf(etas, 1, grid01)
lowercdf <- colMeans(res_$iscontained)
uppercdf <- colMeans(res_$intersects)
# lower / upper 
plot(grid01, lowercdf, type = "l")
lines(grid01, uppercdf)
abline(v = theta_dgp[param], lty = 2)
abline(v = counts[param]/n, lty = 3)

res_2 <- etas_to_lower_upper_cdf_dopar(etas, 1, grid01)
lowercdf2 <- colMeans(res_2$iscontained)
uppercdf2 <- colMeans(res_2$intersects)
# lower / upper 
lines(grid01, lowercdf2, col = "blue", lty = 3)
lines(grid01, uppercdf2, col = "blue", lty = 3)

#### with more categories
# number of observations
n <- 2000
# number of categories
K <- 4
categories <- 1:K
# data 
theta_dgp <- c(0.1, 0.2, 0.3, 0.4)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
counts <- tabulate(X, nbins = K)
cat("Data:", counts, "\n")
#
niterations <- 10000
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
warmup <- 2000
nsubiterations <- 1000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))
#
param <- 4
grid01 <- c(0, seq(from = 0.3, to = 0.5, length.out = 100), 1)
cdfs_ <- etas_to_lower_upper_cdf_dopar(samples_gibbs$etas[subiterations,,], param, grid01)
lowercdf <- colMeans(cdfs_$iscontained)
uppercdf <- colMeans(cdfs_$iscontained)

# lower / upper 
plot(grid01, lowercdf, type = "l")
lines(grid01, uppercdf)
abline(v = theta_dgp[param], lty = 2)
abline(v = counts[param]/n, lty = 3)






