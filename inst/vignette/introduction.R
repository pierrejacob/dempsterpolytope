## This script runs the Gibbs sampler
## on a simulated data, and computes some lower and upper probabilities
## from the output of the Gibbs sampler.

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
counts <- tabulate(X, nbins = K)
print(counts)
## these are the observed data

## run Gibbs sampler to get random polytopes
niterations_gibbs <- 1e3
gibbs_results <- gibbs_sampler(niterations_gibbs, counts)

# the "etas" are stored in gibbs_results$etas, e.g.
# for the etas corresponding to iteration 100
# we get a K x K matrix of etas[k,j] for k,j in 1:K
gibbs_results$etas[100,,]
# from which we can get the vertices of the polytope as
cvxpolytope <- etas_vertices(gibbs_results$etas[100,,])
print(cvxpolytope)
# as we can see this polytope has many vertices, each in the simplex of dimension K

# and we can also get a representation of the polytope 
# as the points x satisfying A x <= b with 
Hrep <- eta_linearconstraints(gibbs_results$etas[100,,])
A <- Hrep$constr
b <- Hrep$rhs
print(A)
print(b)

## from these random polytopes we can estimate lower and upper probabilities 
## associated with assertions of interest.
## Consider the assertion: theta_1 is in the interval [0.2, 0.3].
## To assess whether the sets intersect or are contained in that assertion
## we compute the minimal and maximal value of theta_1 over each set
library(doParallel)
registerDoParallel(cores = detectCores()-2)
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
## for the assertion "theta_1 is in (0.2, 0.3)"
cat(mean(contained_), mean(intersect_), "\n")

## that was for a specific "assertion" on theta_1
## we can also plot the lower and upper CDF associated with any parameter theta_k
## as follows, for component 1.

uppercdf_component1  <- ecdf(min1st)
lowerecdf_component1 <- ecdf(max1st)

# plot lower and upper CDFs for component 1
grid01 <- seq(from = 0, to = 1, length.out = 500)
plot(x = grid01, y = lowerecdf_component1(grid01), type = "l", xlab = expression(theta[1]), ylab = "CDF")
lines(x = grid01, y = uppercdf_component1(grid01), lty = 2)
abline(v = counts[1]/n, lty = 3) # add vertical dotted line at the MLE


