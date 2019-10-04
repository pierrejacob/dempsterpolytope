### This scripts illustrates that inference obtained with or without empty categories are not the same.

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
n <- 10
# number of categories
K <- 3
categories <- 1:K
# data 
theta_dgp <- c(0.5, 0.2, 0.3)
# theta_dgp <- c(0.1, 0.2, 0.3, 0.4)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
freqX <- tabulate(X, nbins = K)
cat("Data:", freqX, "\n")
#
niterations <- 25000
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
samples_gibbs_extended <- gibbs_sampler(niterations, freqX = c(freqX, 0))
names(samples_gibbs)
dim(samples_gibbs$etas_chain)
#
warmup <- 500
nsubiterations <- 5000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))

#
param <- 1
ngrid01 <- 75
grid01 <- seq(from = 0, to = 1, length.out = ngrid01)
etas <- samples_gibbs$etas_chain[subiterations,,]
res_ <- etas_to_lower_upper_cdf_dopar(etas, 1, grid01)
lowercdf <- colMeans(res_$iscontained)
uppercdf <- colMeans(res_$intersects)
# lower / upper 
plot(grid01, lowercdf, type = "l")
lines(grid01, uppercdf)
abline(v = theta_dgp[param], lty = 2)
abline(v = freqX[param]/n, lty = 3)

etas_extended <- samples_gibbs_extended$etas_chain[subiterations,,]
res_extended <- etas_to_lower_upper_cdf_dopar(etas_extended, 1, grid01)
lowercdf_extended <- colMeans(res_extended$iscontained)
uppercdf_extended <- colMeans(res_extended$intersects)
lines(grid01, lowercdf_extended, lty = 2)
lines(grid01, uppercdf_extended, lty = 2)

## add empty category
length(samples_gibbs$Achain)
dim(samples_gibbs$Achain[[1]])

newA <- list()
for (category in 1:K){
  newA[[category]] <- array(NA, dim = c(niterations, freqX[category], K+1))
  for (iter in 1:niterations){
    for (iA in 1:freqX[category]){
      oldA <- samples_gibbs$Achain[[category]][iter,iA,]
      s <- rgamma(1, K, 1)
      w <- rexp(1, 1)
      newA[[category]][iter,iA,] <- c(s * oldA / (s + w), w / (s + w))
    } 
  }
}

length(newA)
dim(newA[[1]])
newA[[1]][1,,]
# apply(newA[[1]][1,,], 1, sum)

## next compute new etas'

new_etas <- array(NA, dim = c(niterations, K + 1, K + 1))
## etas[k,l] = min_u u_l/u_k
for (iter in 1:niterations){
  new_etas[iter,1:K,1:K] <- samples_gibbs$etas_chain[iter,,] 
  new_etas[iter,K+1,] <- Inf
  for (category in 1:K){
    new_etas[iter,category,K+1] <- min(newA[[category]][iter,,K+1]/newA[[category]][iter,,category]) 
  }
}

new_res <- etas_to_lower_upper_cdf_dopar(new_etas, 1, grid01)
new_lowercdf <- colMeans(new_res$iscontained)
new_uppercdf <- colMeans(new_res$intersects)
lines(grid01, new_lowercdf, lty = 3, col = "red")
lines(grid01, new_uppercdf, lty = 3, col = "red")



# samples_gibbs_extended$etas_chain[1,,]

## next, removal of a category
removed_etas <- array(NA, dim = c(niterations, K, K))
for (iter in 1:niterations){
  removed_etas[iter,,] <- samples_gibbs_extended$etas_chain[iter,1:K,1:K]
}
##

res_removed <- etas_to_lower_upper_cdf_dopar(removed_etas[subiterations,,], 1, grid01)
lowercdf_removed <- colMeans(res_removed$iscontained)
uppercdf_removed <- colMeans(res_removed$intersects)
lines(grid01, lowercdf_removed, lty = 2, col = "blue")
lines(grid01, uppercdf_removed, lty = 2, col = "blue")

