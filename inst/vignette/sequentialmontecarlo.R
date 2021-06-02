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
samples_gibbs <- gibbs_sampler(niterations_gibbs, counts)

## from these random polytopes we can estimate lower and upper probabilities 
## associated with assertions of interest
## consider the assertion: theta_1 is in the interval [0.2, 0.3]
## we can encode the interval as a polytope in the simplex of dimension K
## with the instruction:
intervalcvxp <- interval2polytope(K, 1, c(0.2, 0.3))

## then we can compute whether our Gibbs samples intersect or are contained in the 
## polytope corresponding to our interval of interest
## we have to specify some burn-in perod
burnin <- 100
postburn <- niterations_gibbs - burnin
contained_ <- rep(0, postburn)
intersects_ <- rep(0, postburn)
for (index in ((burnin+1):niterations_gibbs)){
  cvxp <- etas2vertices(samples_gibbs$etas[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index-burnin] <- res_[1]
  intersects_[index-burnin] <- res_[2]
}

## this gives us the following lower and upper probabilities
cat(mean(contained_), mean(intersects_), "\n")

## next we can approximate the same quantities but using sequential Monte Carlo
nparticles <- 2^10
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, verbose = FALSE)
## this gives a number of particles, each of which is a random polytope with a weight
## for instance the 10th particle is
samples_smc$etas_particles[10,,]
## with weight
samples_smc$weights[10]

## from which we can compute lower and upper probabilities for the same interval as
contained_smc <- rep(0, nparticles)
intersects_smc <- rep(0, nparticles)
for (iparticle in 1:nparticles){
  cvxp <- etas2vertices(samples_smc$etas_particles[iparticle,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_smc[iparticle] <- res_[1]
  intersects_smc[iparticle] <- res_[2]
}
## the lower and upper probabilities are then approximated as weighted averages over the particles
cat(sum(samples_smc$weights * contained_smc), sum(samples_smc$weights * intersects_smc) , "\n")

## these Monte Carlo estimates provide consistent estimators of the lower and upper
## probabilities of interest. For the Gibbs sampler, as the number of iterations goes to infinity.
## For the SMC sampler, as the number of particles goes to infinity.
