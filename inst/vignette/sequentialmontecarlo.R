## This script runs both the Gibbs sampler and a sequential Monte Carlo algorithm
## and computes some lower and upper probabilities of some assertion
## to compare the results

rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
set.seed(13)

## number of categories
K <- 3
## number of observations
n <- 20
## data-generating parameter value
theta_dgp <- c(0.2, 0.4, 0.4)
## observations, simulated
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
counts <- tabulate(X, nbins = K)
print(counts)
## these are the observed data

## run Gibbs sampler to get random polytopes
niterations_gibbs <- 2e3
gibbs_results <- gibbs_sampler(niterations_gibbs, counts)

## from these random polytopes we can estimate lower and upper probabilities 
## associated with assertions of interest
## consider the assertion: theta_1 is in the interval [0.2, 0.3]

## To assess whether the sets intersect or are contained in that assertion
## we compute the minimal and maximal value of theta_1 over each set

burnin <- 100
postburn <- niterations_gibbs - burnin
etas <- gibbs_results$etas[(burnin+1):niterations_gibbs,,]

min1st <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  lpsolve_over_eta(etas[ieta,,], c(1, rep(0, K-1)))
}
max1st <- foreach(ieta = 1:(dim(etas)[1]), .combine = c) %dopar% {
  -lpsolve_over_eta(etas[ieta,,], c(-1, rep(0, K-1)))
}

## intersection occurs unless max1st < 0.2 or min1st > 0.3
intersect_ <- !((max1st < 0.2) | (min1st > 0.3))
## inclusion occurs when min1st > 0.2 and max1st < 0.3
contained_ <- ((min1st > 0.2) & (max1st < 0.3))

## this gives us the following lower and upper probabilities
cat(mean(contained_), mean(intersect_), "\n")

## next we can approximate the same quantities but using sequential Monte Carlo
nparticles <- 2^10
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, verbose = FALSE)
## this gives a number of particles, each of which is a random polytope with a weight
## for instance the 10th particle is
samples_smc$etas_particles[10,,]
## with weight
samples_smc$weights[10]

min1st_smc <- foreach(ieta = 1:nparticles, .combine = c) %dopar% {
  lpsolve_over_eta(samples_smc$etas_particles[ieta,,], c(1, rep(0, K-1)))
}
max1st_smc <- foreach(ieta = 1:nparticles, .combine = c) %dopar% {
  -lpsolve_over_eta(samples_smc$etas_particles[ieta,,], c(-1, rep(0, K-1)))
}

## intersection occurs unless max1st < 0.2 or min1st > 0.3
intersects_smc <- !((max1st_smc < 0.2) | (min1st_smc > 0.3))
## inclusion occurs when min1st > 0.2 and max1st < 0.3
contained_smc <- ((min1st_smc > 0.2) & (max1st_smc < 0.3))

## the lower and upper probabilities are then approximated as weighted averages over the particles
cat(sum(samples_smc$weights * contained_smc), sum(samples_smc$weights * intersects_smc) , "\n")

## these Monte Carlo estimates provide consistent estimators of the lower and upper
## probabilities of interest. For the Gibbs sampler, as the number of iterations goes to infinity.
## For the SMC sampler, as the number of particles goes to infinity.
