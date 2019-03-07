library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(2)
rm(list = ls())

# number of observations
n <- 50
# number of categories
K <- 4
categories <- 1:K
# data 
## pcol = probability of Y, the column variable, being equal to one
pcol <- 0.4
## prow = probability of X, the row variable, being equal to one
prow <- 0.7
## 
pX <- c(1-prow, prow)
pY <- c(1-pcol, pcol)
# joint probabilities
table_ <- pX %*% t(pY)
table_
# table_ represents
# p_00 p_01
# p_10 p_11
## where p_ij is the probability of X=i, Y=j under independence assumption

X <- sample(x = categories, size = n, replace = TRUE, prob = table_)
freqX <- tabulate(X, nbins = 4)
freqX / sum(freqX)
matrix(freqX / sum(freqX), nrow = 2)

## independence here means theta1 (topleft) * theta4(bottomright) = theta2(bottomleft) * theta3(topright)
## for a given etas, check whether it is compatible with 
##     log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  <= 0 
## and log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  >= 0

# run SMC sampler
nparticles <- 256
samples_smc <- SMC_sampler(nparticles, X, K)
etas <- samples_smc$etas_particles[1,,]
intersect_smc <- foreach(iparticle = 1:nparticles) %dorng% {
  check_intersection_independence(samples_smc$etas_particles[iparticle,,])  
}
intersect1_smc <- sapply(intersect_smc, function(x) x$intersect1)
intersect2_smc <- sapply(intersect_smc, function(x) x$intersect2)
contained1_smc <- sapply(intersect_smc, function(x) x$contained1)
contained2_smc <- sapply(intersect_smc, function(x) x$contained2)

# lower - upper proba on 'theta_1 theta_4 <= theta_2 theta_3
sum(samples_smc$weights * contained1_smc) 
sum(samples_smc$weights * intersect1_smc) 

# lower - upper proba on 'theta_1 theta_4 >= theta_2 theta_3
sum(samples_smc$weights * contained2_smc)
sum(samples_smc$weights * intersect2_smc)

# upper proba on independence 
sum(samples_smc$weights * (intersect1_smc & intersect2_smc))
