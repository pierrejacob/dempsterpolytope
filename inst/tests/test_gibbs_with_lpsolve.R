library(montecarlodsm)
set.seed(1)
rm(list = ls())
# install.packages("lpSolveAPI")
library(lpSolveAPI)

#
## Let's generate some data 
## number of categories
K <- 5
## number of observations
n <- 50
theta_dgp <- rexp(K)
theta_dgp <- theta_dgp / sum(theta_dgp)
## observations
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K) + rep(1, K)
print(freqX)
## (note that each count was incremented by one to ensure non zero entries, for simplicity)



niterations_gibbs <- 200

## run Gibbs sampler with graph tools
set.seed(1)
samples_gibbs <- gibbs_sampler_graph(niterations_gibbs, freqX)
etas <- samples_gibbs$etas_chain[niterations_gibbs,,]
etas
## run Gibbs sampler with lpSolve
set.seed(1)
samples_gibbs_lpsolve <- gibbs_sampler_lp(niterations_gibbs, freqX)
etas_lpsolve <- samples_gibbs_lpsolve$etas_chain[niterations_gibbs,,]
etas_lpsolve
#
summary(as.numeric(abs(etas - etas_lpsolve)))

## 
library(microbenchmark)
microbenchmark(samples_gibbs_graph <- gibbs_sampler_graph(niterations_gibbs, freqX),
               samples_gibbs_lp <- gibbs_sampler_lp(niterations_gibbs, freqX),
               times = 10)
# 


