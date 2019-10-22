### This scripts illustrates that inference obtained with or without empty categories are not the same.
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of categories
K <- 3
categories <- 1:K
# data 
freqX <- c(4,3,2)
cat("Data:", freqX, "\n")
##
niterations <- 10000
samples_gibbs_K3 <- gibbs_sampler(niterations, freqX)
samples_gibbs_K4 <- gibbs_sampler(niterations, c(freqX, 0))
##
## burnin and thinning
warmup <- 1000
nsubiterations <- 2000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))
etas_K3 <- samples_gibbs_K3$etas_chain[subiterations,,]
etas_K4 <- samples_gibbs_K4$etas_chain[subiterations,,]

## compare empirical lower/upper CDF of first parameter
ngrid01 <- 100
grid01 <- seq(from = 0, to = 1, length.out = ngrid01)

res_K31 <- etas_to_lower_upper_cdf_dopar(etas_K3, 1, grid01)
lowercdf_K31 <- colMeans(res_K31$iscontained)
uppercdf_K31 <- colMeans(res_K31$intersects)

res_K41 <- etas_to_lower_upper_cdf_dopar(etas_K4, 1, grid01)
lowercdf_K41 <- colMeans(res_K41$iscontained)
uppercdf_K41 <- colMeans(res_K41$intersects)

# lower / upper 
g1 <- qplot(x = grid01, y = lowercdf_K31, linetype = "3", geom = "line") + xlab(expression(theta[1])) + ylab("CDF") +
  geom_line(aes(y = uppercdf_K31, linetype = "3")) + 
  geom_line(aes(y = lowercdf_K41, linetype = "4")) + 
  geom_line(aes(y = uppercdf_K41, linetype = "4")) + scale_linetype(name = "# categories: ")
g1

ggsave(filename = "inst/reproduce/emptycategory.theta1.pdf", plot = g1, width = 7, height = 5)
###
res_K32 <- etas_to_lower_upper_cdf_dopar(etas_K3, 2, grid01)
lowercdf_K32 <- colMeans(res_K32$iscontained)
uppercdf_K32 <- colMeans(res_K32$intersects)

res_K42 <- etas_to_lower_upper_cdf_dopar(etas_K4, 2, grid01)
lowercdf_K42 <- colMeans(res_K42$iscontained)
uppercdf_K42 <- colMeans(res_K42$intersects)

# lower / upper 
g2 <- qplot(x = grid01, y = lowercdf_K32, linetype = "3", geom = "line") + xlab(expression(theta[2])) + ylab("CDF") +
  geom_line(aes(y = uppercdf_K32, linetype = "3")) + 
  geom_line(aes(y = lowercdf_K42, linetype = "4")) + 
  geom_line(aes(y = uppercdf_K42, linetype = "4")) + scale_linetype(name = "# categories: ")
g2

ggsave(filename = "inst/reproduce/emptycategory.theta2.pdf", plot = g2, width = 7, height = 5)



### Next we see if we can retrieve the same plot by 
## taking the samples obtained with K = 3 and extending them to K = 4
## taking the samples obtained with K = 4 and reducing them to K = 4

## first, extend from K = 3
pts_extended_from_K3 <- list()
for (category in 1:K){
  pts_extended_from_K3[[category]] <- array(NA, dim = c(niterations, freqX[category], K+1))
  for (iter in 1:niterations){
    for (iA in 1:freqX[category]){
      oldA <- samples_gibbs_K3$Achain[[category]][iter,iA,]
      s <- rgamma(1, K, 1)
      w <- rexp(1, 1)
      pts_extended_from_K3[[category]][iter,iA,] <- c(s * oldA / (s + w), w / (s + w))
    } 
  }
}
## next compute new etas'
etas_K4_from_K3 <- array(NA, dim = c(niterations, K + 1, K + 1))
## etas[k,l] = min_u u_l/u_k
for (iter in 1:niterations){
  etas_K4_from_K3[iter,1:K,1:K] <- samples_gibbs_K3$etas_chain[iter,,] 
  etas_K4_from_K3[iter,K+1,] <- Inf
  for (category in 1:K){
    etas_K4_from_K3[iter,category,K+1] <- min(pts_extended_from_K3[[category]][iter,,K+1]/pts_extended_from_K3[[category]][iter,,category]) 
  }
}


## next, removal of a category in K = 4
## that's simpler...
etas_K3_from_K4 <- array(NA, dim = c(niterations, K, K))
for (iter in 1:niterations){
  etas_K3_from_K4[iter,,] <- samples_gibbs_K4$etas_chain[iter,1:K,1:K]
}
##

alternative_res_K3 <- etas_to_lower_upper_cdf_dopar(etas_K3_from_K4, 1, grid01)
alternative_lowercdf_K3 <- colMeans(alternative_res_K3$iscontained)
alternative_uppercdf_K3 <- colMeans(alternative_res_K3$intersects)
alternative_res_K4 <- etas_to_lower_upper_cdf_dopar(etas_K4_from_K3, 1, grid01)
alternative_lowercdf_K4 <- colMeans(alternative_res_K4$iscontained)
alternative_uppercdf_K4 <- colMeans(alternative_res_K4$intersects)

# lower / upper 
galt <- qplot(x = grid01, y = alternative_lowercdf_K3, linetype = "3", geom = "line") + xlab(expression(theta[1])) + ylab("CDF") +
  geom_line(aes(y = alternative_uppercdf_K3, linetype = "3")) + 
  geom_line(aes(y = alternative_lowercdf_K4, linetype = "4")) + 
  geom_line(aes(y = alternative_uppercdf_K4, linetype = "4")) + scale_linetype(name = "# categories: ")
galt
### this should look like g1
