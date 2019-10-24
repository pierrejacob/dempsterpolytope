## This scripts runs the Gibbs sampler 
## on counts in 3 categories
## and shows resulting polytopes in plots

library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

## number of observations
n <- 50
## number of categories
K <- 3
categories <- 1:K
## data-generating parameter value
theta_dgp <- c(0.2, 0.4, 0.4)
## observations, simulated
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_dgp)
## frequencies
counts <- tabulate(X, nbins = K)
print(counts/sum(counts))
##
## encompassing triangle has three vertices
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")

###
niterations <- 1000
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts, theta_0 = c(0.8,0.1,0.1))
(proc.time() - pct)[3]
##
etas <- samples_gibbs$etas[501:1000,,]

## and overlay them in plot
df.polytope <- data.frame()
for (iteration in 1:200){
  etas_iteration <- etas[iteration,,]
  etascvxp <- etas2cvxpolytope(etas_iteration)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                               iteration = iteration))
}
g <- ggplot_triangle(v_cartesian) +
  geom_polygon(data=df.polytope %>% filter(iteration >= 100), aes(x = x, y = y, group = iteration), alpha = .3)
g

###
## Suppose we want to make assertion on "X_{n+1} = k"
k <- 2
## random sets associated with the assertion
## random point in the simplex
pointv <- rexp(K, 1)
pointv <- pointv / sum(pointv)
pointv_cart <- barycentric2cartesian(pointv, v_cartesian)
## region Sigma(v, k) = {theta in simplex: theta_{ell} / theta_k <= v_ell / v_k}
## the constraints are on the first K-1 coordinates
## the first ones say that the feasible set is within the simplex
predictive_region <- function(pointv, k){
  K_ <- length(pointv)
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K_-1, K_-1))
  b <- c(1, rep(0, K_-1))
  # then the extra constraints come from v: the points x should satisfy x_ell / x_k <= v_ell / v_k for all ell not equal to k 
  for (j in setdiff(categories, k)){
    ccc <- rep(0, K_)
    ccc[k] <- -pointv[j]/pointv[k]
    ccc[j] <- 1
    cc <- ccc - ccc[K_]
    b <- c(b, -ccc[K_])
    A <- rbind(A, matrix(cc[1:(K_-1)], nrow = 1))
  }
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  return(list(constr = constr))
}

pointv_constr <- predictive_region(pointv, k)
hrepr <- rcdd::makeH(pointv_constr$constr$constr, pointv_constr$constr$rhs)
vrepr <- rcdd::q2d(rcdd::scdd(rcdd::d2q(hrepr))$output)
vrepr <- cbind(vrepr[,-c(1,2)], 1- apply(vrepr[,-c(1,2)], 1, sum))
vrepr_cart <- t(apply(vrepr, 1, function(x) barycentric2cartesian(x, v_cartesian)))
## order polygon
average_ <- colMeans(vrepr_cart)
o_ <- order(apply(sweep(vrepr_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
vrepr_cart <- vrepr_cart[o_,]


## compute intersections / containment in feasible sets
whichcontained  <- rep(FALSE, dim(etas)[1])
whichintersects <- rep(FALSE, dim(etas)[1])
for (iteration in 1:dim(etas)[1]){
  cvx_gibbs <- etas2cvxpolytope(etas[iteration,,])
  res_ <- compare_polytopes(cvx_gibbs, pointv_constr)
  whichcontained[iteration] <- res_[1]
  whichintersects[iteration] <- res_[2]
}

## color differently polytopes which intersects, are contained, or neither
df.polytope <- data.frame()
for (iteration in 1:dim(etas)[1]){
  etas_iteration <- etas[iteration,,]
  etascvxp <- etas2cvxpolytope(etas_iteration)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  type <- whichcontained[iteration] + whichintersects[iteration]
  df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                               iteration = iteration, type = type))
}
g <- ggplot_triangle(v_cartesian) +
  geom_polygon(data=df.polytope %>% filter(iteration >= 100), aes(x = x, y = y, group = iteration, fill = factor(type)), alpha = .3)
g + geom_polygon(data=data.frame(x = vrepr_cart[,1], y = vrepr_cart[,2]), fill = "black", alpha = 0.2) + 
  geom_point(data = data.frame(x = pointv_cart[1], y = pointv_cart[2]), col = "red")

## Now do the computation for many points v in the simplex
nv <- 1e1
lower_prob <- 0
upper_prob <- 0
for (iteration in 1:dim(etas)[1]){
  pointvs <- matrix(rexp(K*nv, 1), nrow = nv)
  pointvs <- t(apply(pointvs, 1, function(v) v / sum(v)))
  whichcontained <- rep(0, nv)
  whichintersects <- rep(0, nv)
  cvx_gibbs <- etas2cvxpolytope(etas[iteration,,])
  for (iv in 1:nv){
    pointv_constr <- predictive_region(pointvs[iv,], k)
    res_ <- compare_polytopes(cvx_gibbs, pointv_constr)
    whichcontained[iv] <- res_[1]
    whichintersects[iv] <- res_[2]
  }
  ## lower probability 
  lower_prob <- lower_prob + mean(whichcontained)/dim(etas)[1]
  ## upper probability 
  upper_prob <- upper_prob + mean(whichintersects)/dim(etas)[1]  
}

cat("assertion: next observation will be in category", k, "\n")
cat("probabilities", lower_prob, upper_prob, "\n")
cat("compared to empirical frequencies:\n")
print(counts/sum(counts))

