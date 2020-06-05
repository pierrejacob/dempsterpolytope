## This scripts tests the method to add/remove empty categories 
## i.e. how to modify the 'uniform points' denoted by u_n in the paper

rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
set.seed(4)

niterations <- 25000
counts_K4 <- c(7,5,8,0)
samples_gibbs_K4 <- gibbs_sampler(niterations = niterations, counts = counts_K4, theta_0 = c(0.7,0.1,0.1,0.1))
#
counts_K3 <- c(7,5,8)
samples_gibbs_K3 <- gibbs_sampler(niterations = niterations, counts = counts_K3, theta_0 = c(0.8,0.1,0.1))


## look at a_n where n is first index in I_1
# hist(samples_gibbs_K3$Us[[1]][,1,1], nclass = 100)
# hist(samples_gibbs_K4$Us[[1]][,1,1], nclass = 100)

## if we add an empty category to results obtained with K = 3...
## for pts in 1st category 
dim(samples_gibbs_K3$Us[[1]])
new_Us1 <- array(NA, dim = c(dim(samples_gibbs_K3$Us[[1]])[1], counts_K3[1], 4))
for (n in 1:counts_K3[1]){
  pts_k3 <- samples_gibbs_K3$Us[[1]][,n,]
  pts_k3_augmented <- t(apply(pts_k3, 1, function(v){
    sum_ <- rgamma(1, 3, 1)
    y <- c(v * sum_, rexp(1, 1))
    return(y / sum(y))
  } ))
  new_Us1[,n,] <- pts_k3_augmented
}
hist(new_Us1[100:niterations,2,2], prob = TRUE, nclass = 50)
hist(samples_gibbs_K4$Us[[1]][100:niterations,2,2], prob = TRUE, nclass = 50, add = TRUE, col = rgb(1,0,0,0.5))

## if we remove the empty category to results obtained with K = 4...

dim(samples_gibbs_K4$Us[[1]])
new_Us1 <- array(NA, dim = c(dim(samples_gibbs_K3$Us[[1]])[1], counts_K4[1], 3))
for (n in 1:counts_K4[1]){
  pts_k4 <- samples_gibbs_K4$Us[[1]][,n,]
  pts_k4_reduced <- t(apply(pts_k4, 1, function(v){
    y <- v[-4]
    return(y / sum(y))
  } ))
  new_Us1[,n,] <- pts_k4_reduced
}
hist(new_Us1[100:niterations,3,3], prob = TRUE, nclass = 50)
hist(samples_gibbs_K3$Us[[1]][100:niterations,3,3], prob = TRUE, nclass = 50, add = TRUE, col = rgb(1,0,0,0.5))
