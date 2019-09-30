library(montecarlodsm)
set.seed(1)
rm(list = ls())

## Let's generate some data 
## number of categories
K <- 5
## number of observations
n <- 50
## data-generating value
theta_dgp <- c(0.2, 0.4, 0.2, 0.1, 0.1)
## observations
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K)
print(freqX)

## run Gibbs sampler to get random polytopes
niterations_gibbs <- 1e3
samples_gibbs <- gibbs_sampler(niterations_gibbs, freqX)

# the "etas" are stored in samples_gibbs$etas_chain, e.g.
# for the etas corresponding to iteration 100
# we get a K x K matrix of etas[k,j] for k,j in 1:K
samples_gibbs$etas_chain[100,,]
# from which we can get the vertices of the polytope as
cvxpolytope <- etas2cvxpolytope(samples_gibbs$etas_chain[100,,])
print(cvxpolytope$vertices_barcoord)
# as we can see this polytope has 10 vertices, each in the simplex of dimension 4

# and we can also get a representation of the polytope 
# as the points x satisfying A x <= b with 
A <- cvxpolytope$constr$constr
b <- cvxpolytope$constr$rhs
print(A)
print(b)

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
# library(doParallel)
# registerDoParallel(cores = detectCores()-2)
# rescvx <- foreach(index = ((burnin+1):niterations_gibbs), .combine = rbind) %dopar% {
#   cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
#   res_ <- compare_polytopes(cvxp, intervalcvxp)
#   res_
# }
# contained_ <- rescvx[,1]
# intersects_ <- rescvx[,2]
for (index in ((burnin+1):niterations_gibbs)){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index-burnin] <- res_[1]
  intersects_[index-burnin] <- res_[2]
}

## this gives us the following lower and upper probabilities
cat(mean(contained_), mean(intersects_), "\n")

## next we can approximate the same quantities but using a different Monte Carlo technique
## called "sequential Monte Carlo"
nparticles <- 2^10
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, verbose = FALSE)
## this gives a number of particle, each of which is a random polytope with a weight
## for instance the 10th particle is
samples_smc$etas_particles[10,,]
## with weight
samples_smc$weights[10]

## from which we can compute lower and upper probabilities for the same interval as
contained_smc <- rep(0, nparticles)
intersects_smc <- rep(0, nparticles)
for (iparticle in 1:nparticles){
  cvxp <- etas2cvxpolytope(samples_smc$etas_particles[iparticle,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_smc[iparticle] <- res_[1]
  intersects_smc[iparticle] <- res_[2]
}
## the lower and upper probabilities are then approximated as weighted averages over the particles
cat(sum(samples_smc$weights * contained_smc), sum(samples_smc$weights * intersects_smc) , "\n")

## these Monte Carlo estimates provide consistent estimators of the lower and upper
## probabilities of interest. For the Gibbs sampler, as the number of iterations goes to infinity.
## For the SMC sampler, as the number of particles goes to infinity.
