## This script runs the Gibbs sampler
## on a simulated data, and computes some lower and upper probabilities
## from the output of the Gibbs sampler

library(dempsterpolytope)
set.seed(1)
rm(list = ls())

## number of categories
K <- 5
## number of observations
n <- 50
## data-generating parameter value
theta_dgp <- c(0.2, 0.4, 0.2, 0.1, 0.1)
## observations, simulated
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
freqX <- tabulate(X, nbins = K)
print(freqX)
## these are the observed data

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
# as we can see this polytope has many vertices, each in the simplex of dimension K

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

## then we can compute whether our polytopes intersect or are contained in the 
## polytope corresponding to our interval of interest
## we have to specify some burn-in perod
burnin <- 100
postburn <- niterations_gibbs - burnin
contained_ <- rep(0, postburn)
intersects_ <- rep(0, postburn)
for (index in ((burnin+1):niterations_gibbs)){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index-burnin] <- res_[1]
  intersects_[index-burnin] <- res_[2]
}

## this gives us the following lower and upper probabilities
## for the assertion "theta_1 is in (0.2, 0.3)"
cat(mean(contained_), mean(intersects_), "\n")

## that was for a specific "assertion" on theta_1
## we can also plot the lower and upper CDF associated with any parameter theta_k
## as follows

## we are going to speed things up a little bit with parallel computation
library(doParallel)
registerDoParallel(cores = detectCores()-2)

grid01 <- seq(from = 0, to = 1, length.out = 50)
loweruppercdf <- etas_to_lower_upper_cdf_dopar(samples_gibbs$etas_chain[(burnin+1):niterations_gibbs,,], 1, grid01)
lowercdf_ <- colMeans(loweruppercdf$iscontained)
uppercdf_ <- colMeans(loweruppercdf$intersects)

# plot lower and upper CDFs
plot(x = grid01, y = lowercdf_, type = "l", xlab = expression(theta[1]), ylab = "CDF")
lines(x = grid01, y = uppercdf_)
abline(v = freqX[1]/n, lty = 3) # add vertical dotted line at the MLE