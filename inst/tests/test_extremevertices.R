## This implements Art's idea
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(2)
##

K <- 4

counts <- c(3,2,1,5)
N <- sum(counts)

niterations <- 5e3
gibbs_results <- dempsterpolytope::gibbs_sampler(niterations, counts)
names(gibbs_results)


burnin <- 1e3
etas <- gibbs_results$etas[(burnin+1):niterations,,]
## vertices with maximum component 'param'
objvec <- rep(0, K)
param <- 4
objvec[param] <- 1
max_param <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  -lpsolve_over_eta(etas[ieta,,], -objvec)
}

# sample from Dirichlet distribution supposed to match
dirisamples <- t(gtools::rdirichlet(1e4, alpha = counts + objvec))
#
hist(max_param, nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[param,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
