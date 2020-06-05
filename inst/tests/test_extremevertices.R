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
etas_iteration <- gibbs_results$etas[burnin,,]
etascvxp <- etas2cvxpolytope(etas_iteration)
# ?etas2cvxpolytope

## vertex with maximum 1st component
max(etascvxp$vertices_barcoord[,1])

maxv1 <- sapply((burnin+1):niterations, function(index){
  polytope <- etas2cvxpolytope(gibbs_results$etas[index,,])$vertices_barcoord
  polytope[which.max(polytope[,1]),]
  })

# sample from Dirichlet distribution supposed to match
dirisamples <- t(gtools::rdirichlet(1e4, alpha = counts + c(1,0,0,0)))
#
par(mfrow = c(2,2))
hist(maxv1[1,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[1,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(maxv1[2,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[2,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(maxv1[3,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[3,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(maxv1[4,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[4,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)


minv1 <- sapply((burnin+1):niterations, function(index){
  polytope <- etas2cvxpolytope(gibbs_results$etas[index,,])$vertices_barcoord
  polytope[which.min(polytope[,1]),]
})

dirisamples <- t(gtools::rdirichlet(1e4, alpha = counts + c(0,1,1,1)))
#
par(mfrow = c(2,2))
hist(minv1[1,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[1,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(minv1[2,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[2,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(minv1[3,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[3,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)
hist(minv1[4,], nclass=30, xlim=c(0,1), prob=TRUE)
hist(dirisamples[4,], add=TRUE, prob=TRUE, col = rgb(1,0,.5,0.5), nclass=30)


