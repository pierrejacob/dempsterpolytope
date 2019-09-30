library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations
n <- 22
# number of categories
K <- 4
categories <- 1:K
# data 
# freqX <- c(70,150,0)
# theta_dgp <- c(0.3, 0.3, 0.4)
# X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)

X <- c(categories, sample(x = categories, size = n - K, replace = TRUE))
freqX <- tabulate(X)
freqX

### SMC sampler
nparticles <- 2^5
smc_res <- SMC_sampler(nparticles, X, K)
sum(smc_res$normcst)
resamplingtimes <- smc_res$resamplingtimes

# smc_res <- SMC_sampler(nparticles, X, K, resamplingtimes = resamplingtimes)

NREP <- 100
nparticles <- 2^5
smc_result <- SMC_sampler(nparticles, X, K)

smcsampler_results <- foreach(irep = 1:NREP) %dopar% {
  SMC_sampler(nparticles, X, K, resamplingtimes = smc_result$resamplingtimes)
}
normcsts <- sapply(smcsampler_results, function(x) sum(x$normcst))

var(normcsts)
hist(normcsts)

## from there, deduce what the meeting time would be like
meanaccepts <- foreach (i = 1:NREP, .combine = c) %dopar% {
  Zstart <- normcsts[i]
  othercsts <- normcsts[-i]
  mean(pmin(1, exp(othercsts - Zstart)))
}
# sample from mixture of Geometric
fake_meetings <- 1 + rgeom(length(meanaccepts), prob = meanaccepts)
mean(fake_meetings > 2)
summary(fake_meetings)
hist(fake_meetings)

# ## check against real meeting times 
cchains <- foreach(irep = 1:NREP) %dorng% {
  coupled_chains(nparticles, X, K, resamplingtimes)
}
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
summary(meetingtimes)
mean(meetingtimes > 2)
hist(meetingtimes)
## seems like an accurate approximation!

# from there we can take k = 2, m = 5 for instance

h <- function(etas) unlist(check_intersection_independence(etas))
h(smc_result$etas_particles[1,,])
# coupled_chains(nparticles, X, K, resamplingtimes, k = 2, m = 5, h = h)

#### now compare results for independence assumption
#### run standard SMC sampler
nparticles <- 2^8
pct <- proc.time()
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.5)
print(proc.time()-pct)
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

## these sum to one
# sum(samples_smc$weights * contained1_smc)  + sum(samples_smc$weights * contained2_smc)  + sum(samples_smc$weights * (intersect1_smc & intersect2_smc))


## Now compared to unbiased estimators
nparticles <- 32
preliminary_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75)
h <- function(etas) unlist(check_intersection_independence(etas))
NREP <- 10
cchains <- foreach(irep = 1:NREP) %dorng% {
  coupled_chains(nparticles, X, K, resamplingtimes = preliminary_smc$resamplingtimes, k = 1, m = 2, h = h)
}
# 
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
summary(meetingtimes)
sapply(cchains, function(x) x$iteration)
uestimators <- t(sapply(cchains, function(x) x$uestimator))
colMeans(uestimators)
1.96 * apply(uestimators, 2, sd)
