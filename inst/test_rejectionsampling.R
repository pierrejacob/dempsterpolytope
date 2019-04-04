library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations 
# (don't set it to be too large, because rejection sampler would become very slow; values of <= 5,6 are OK)
n <- 6
# number of categories
K <- 3
categories <- 1:K
# data 
# freqX <- c(70,150,0)
# theta_dgp <- c(0.3, 0.3, 0.4)
# X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
X <- c(categories, sample(x = categories, size = n - K, replace = TRUE))
freqX <- tabulate(X)
freqX
##
## encompassing triangle has three vertices
if (K == 3){
  v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
  cols <- c("red", "green", "blue")
}

###
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
##
samples_gibbs$etas_chain[100,,]

### Rejection sampling
rejectionsampler(X, K)

# or manually, for K = 3
# check_cst_manual <- function(etas){
#   satis <- TRUE
#   for (k in 1:K){
#     notk <- setdiff(1:K, k)
#     for (j in notk){
#       l <- setdiff(1:K, c(k,j))
#       if (etas[k,j] < 1/(etas[j,k])){
#         satis <- FALSE
#       }
#       if (etas[k,j] < 1/(etas[j,l] * etas[l,k])){
#         satis <- FALSE
#       }
#     }
#   }
#   return(satis)
# }

NREP <- 1e3
etas_ <- foreach(irep = 1:NREP) %dorng% {
  rejectionsampler(X, K)
}

## to estimate the volume of accept region
## we can use the fact that the number of attempts before a success
## is geometric 
hist(sapply(etas_, function(x) x$nattempts-1))
## so an unbiased estimator of 1/Z is 
mean(sapply(etas_, function(x) x$nattempts))
## and an unbiased estimator of Z is 
1/(NREP/(NREP-1) * mean(sapply(etas_, function(x) x$nattempts-1)) + 1)
## the MLE is 
1/(mean(sapply(etas_, function(x) x$nattempts-1)) + 1)


## check with SMC
nparticles <- 2^12
smc_res <- SMC_sampler(nparticles, X, K)
## log normalizing constant
exp(sum(smc_res$normcst))

## seems to agree

# etas_accepted <- etas_[which(sapply(etas_, function(x) x$sat))]
# if (K == 3){
#   rejection.df.polytope <- data.frame()
#   for (iteration in 1:length(etas_accepted)){
#     etas <- etas_accepted[[iteration]]$etas
#     etascvxp <- etas2cvxpolytope(etas)
#     ## convert coordinates to cartesian
#     vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
#     # order vertices according to angles
#     average_ <- colMeans(vertices_cart)
#     o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
#     vertices_cart <- vertices_cart[o_,]
#     rejection.df.polytope <- rbind(rejection.df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
#                                                                      iteration = iteration))
#   }
#   g <- ggplot_triangle(v_cartesian) +
#     geom_polygon(data=rejection.df.polytope, aes(x = x, y = y, group = iteration), alpha = .3)
#   g
# }  
# # ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.polytopeS.pdf", plot = g, width = 7, height = 7)


##
niterations <- 1100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
elapsed <- (proc.time() - pct)[3]
elapsed
burnin <- 100

param <- 1
interval <- c(0.1, 0.2)
intervalcvxp <- interval2polytope(K, param, interval)
# 
postburn <- niterations - burnin
contained_ <- rep(0, postburn)
intersects_ <- rep(0, postburn)
for (index in ((burnin+1):niterations)){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index-burnin] <- res_[1]
  intersects_[index-burnin] <- res_[2]
}

# lower/upper probabilities, post burn-in
cat(mean(contained_), mean(intersects_), "\n")
# for the interval
interval
# on parameter
param
# equal to 
(freqX/sum(freqX))[param]

## now with rejection sampling
nrs <- 1e3
etas_rs <- foreach(irep = 1:nrs) %dorng% {
  rejectionsampler(X,K)
}
contained_rs <- rep(0, nrs)
intersects_rs <- rep(0, nrs)
for (index in 1:nrs){
  cvxp <- etas2cvxpolytope(etas_rs[[index]]$etas)
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_rs[index] <- res_[1]
  intersects_rs[index] <- res_[2]
}

# lower/upper probabilities, obtained with rejection sampling
cat(mean(contained_rs), mean(intersects_rs), "\n")
# versus the Gibbs output 
cat(mean(contained_), mean(intersects_), "\n")

# seems to agree
