library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations
n <- 20
# number of categories
K <- 4
categories <- 1:K
# data 
pX <- c(0.2, 0.8)
pY <- c(0.4, 0.6)
table_ <- pX %*% t(pY)
table_[1,1] * table_[2,2]
table_[2,1] * table_[1,2]
table_
X <- sample(x = categories, size = n, replace = TRUE, prob = table_)
# X <- sample(x = c(rep(1, 8), rep(2, 32), rep(3, 12), rep(4, 48)), size = n, replace = FALSE)
freqX <- tabulate(X, nbins = 4)
freqX / sum(freqX)
matrix(freqX / sum(freqX), nrow = 2)
## independence means theta1 (topleft) * theta4(bottomright) = theta2(bottomleft) * theta3(topright)

#
niterations <- 2000
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
etas <- samples_gibbs$etas_chain[81,,]

# now we are interested in the set of theta such that theta_ell / theta_k < etas[k,ell]
# etas2cvxpolytope(etas)$vertices_barcoord
#
## Now suppose that we want to do a change of coordinates,
## t_k = log(theta_k) - (1/K) sum(log(theta_k))
## so that sum_k t_k = 0
## then theta_ell / theta_k < c
## is equivalent to log(theta_ell) - log(theta_k) < log(c)
## i.e. t_ell - t_k < log(c)
## 
# K_ <- K
# A <- matrix(0, nrow = K_*(K_-1), ncol = K_-1)
# b <- rep(0, K*(K_-1))
# # then the extra constraints come from etas
# index <- 1
# for (d in categories){
#   for (j in setdiff(categories, d)){
#     # cccc (wA wB wC ... )' = 0
#     ccc <- rep(0, K_)
#     ccc[d] <- -1
#     ccc[j] <- +1
#     cc <- ccc - ccc[K_]
#     b[index] <- log(etas[d,j])
#     A[index,] <- cc[1:(K_-1)]
#     index <- index + 1
#   }
# }
# constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
# # get vertices of polytope
# vertices_barcoord <- hitandrun::findVertices(constr)
# # then add last coordinate, so that entries sum to zero again
# vertices_barcoord <- cbind(vertices_barcoord, - apply(vertices_barcoord, 1, sum))
# # convert back to original scale
# t(apply(vertices_barcoord, 1, function(v) exp(v)/(sum(exp(v)))))
# etas2cvxpolytope(etas)$vertices_barcoord
# ## agreement on the vertices, order has changed

check_intersection_independence <- function(etas){
  ## check intersection between polypote of theta s.t. theta_j / theta_d < etas[d,j] 
  ## with independence assumption
  ## log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3) = 0
  ## i.e. log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  <= 0
  ## and  log(theta_1) + log(theta_4) - log(theta_2) - log(theta_3)  >= 0
  K_ <- dim(etas)[1]
  A <- matrix(0, nrow = K_*(K_-1), ncol = K_-1)
  b <- rep(0, K*(K_-1))
  # then the extra constraints come from etas
  index <- 1
  for (d in categories){
    for (j in setdiff(categories, d)){
      # cccc (wA wB wC ... )' = 0
      ccc <- rep(0, K_)
      ccc[d] <- -1
      ccc[j] <- +1
      cc <- ccc - ccc[K_]
      b[index] <- log(etas[d,j])
      A[index,] <- cc[1:(K_-1)]
      index <- index + 1
    }
  }
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  Aindep <- matrix(c(1,-1,-1,1), nrow = 1, byrow = T)
  Aindep <- Aindep[,1:(K_-1),drop=F] - Aindep[,K_]
  bindep <- 0
  ncst <- nrow(constr$constr)
  ## test intersection: 
  intersectconstr <- list(constr = rbind(constr$constr, Aindep), rhs = c(constr$rhs, bindep), dir = rep("<=", ncst+1))
  intersect1 <- !inherits(try(hitandrun::findVertices(intersectconstr), silent = T), "try-error")
  intersectconstr <- list(constr = rbind(constr$constr, -Aindep), rhs = c(constr$rhs, bindep), dir = rep("<=", ncst+1))
  intersect2 <- !inherits(try(hitandrun::findVertices(intersectconstr), silent = T), "try-error")
  return(intersect1 && intersect2)
}

intersect_ <- as.numeric(foreach(iteration = 1:niterations, .combine = c) %dorng% {
  check_intersection_independence(samples_gibbs$etas_chain[iteration,,])  
})
mean(intersect_)

##
nparticles <- 256
samples_smc <- SMC_sampler(nparticles, X, K)
intersect_smc <- as.numeric(foreach(iparticle = 1:nparticles, .combine = c) %dorng% {
  check_intersection_independence(samples_smc$etas_particles[iparticle,,])  
})
sum(samples_smc$weights * intersect_smc)



