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
matrix(freqX, nrow = 2)

## independence here means theta1 (topleft) * theta4(bottomright) = theta2(bottomleft) * theta3(topright)
## for a given etas, check whether it is compatible with 
##     log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  <= 0 
## and log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  >= 0

# run SMC sampler
nparticles <- 128
h <- function(etas) unlist(check_intersection_independence(etas))
pct <- proc.time()
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, h = h)
print((proc.time() - pct)[3])

p <- samples_smc$hestimator[2]
q <- samples_smc$hestimator[4]
r <- 1-p-q
cat(p, q, r, "\n")

## unbiased estimators
NREP <- 20
pct <- proc.time()
nparticles <- 2^7
smc_result <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, h = h)
smcsampler_results <- foreach(irep = 1:NREP) %dorng% {
  SMC_sampler(nparticles, X, K, resamplingtimes = smc_result$resamplingtimes, h = h)
}
normcsts <- sapply(smcsampler_results, function(x) sum(x$normcst))

## deduce resulting meeting times
meanaccepts <- foreach (i = 1:NREP, .combine = c) %dopar% {
  Zstart <- normcsts[i]
  othercsts <- normcsts[-i]
  mean(pmin(1, exp(othercsts - Zstart)))
}
# sample from mixture of Geometric
fake_meetings <- 1 + rgeom(length(meanaccepts), prob = mean(meanaccepts))
k <- 2
m <- 2*k

cchains <- foreach(irep = 1:NREP) %dorng% {
  coupled_chains(nparticles, X, K, resamplingtimes = smc_result$resamplingtimes, k = k, m = m, h = h)
}

uestimators <- t(sapply(cchains, function(x) x$uestimator))
colMeans(uestimators) - 1.96 * apply(uestimators, 2, sd) / sqrt(NREP)
colMeans(uestimators) + 1.96 * apply(uestimators, 2, sd) / sqrt(NREP)
print((proc.time() - pct)[3])

hist(sapply(cchains, function(x) x$meetingtime))

#