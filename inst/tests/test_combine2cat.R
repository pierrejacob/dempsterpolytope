library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(2)
rm(list = ls())
##
v_cartesian <- list(c(0,0), c(1,0), c(1/2, sin(pi/3)))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
gtriangle <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
gtriangle <- gtriangle + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
K <- 3

## Data
counts <- c(2,1,3)
N <- sum(counts)
## rejection sampler: rejectionsampler(counts)
## Gibbs sampler: samples_gibbs <- gibbs_sampler(niterations = 100, counts = counts)

## given two counts, returns interval of feasible set 
## interval is returned as pair of values in (0,1), each measuring some length from the left endpoint
## i.e. the 'theta2' coordinate of a point (theta1, theta2) in barycentric coordinate 
sample_feasible_set <- function(count1, count2){
  gamma1 <- rgamma(n = 1, shape = count1, rate = 1)
  gamma2 <- rgamma(n = 1, shape = count2, rate = 1)
  wi <- rexp(1, 1)
  interval_ <- c(gamma1/(gamma1+wi+gamma2), (gamma1+wi)/(gamma1+wi+gamma2))
  return(interval_)
}
## transforms interval in 'theta2' coordinate into etas_1->2 and etas_{2->1}
## returns c(etas_1->2, etas_2->1)
interval_to_etas <- function(interval){
  return(c(interval[2] / (1 - interval[2]), (1 - interval[1]) / interval[1]))
}
sample_feasible_set(100, 2)
sample_feasible_set(2, 100)

interval_to_etas(sample_feasible_set(counts[1], counts[2]))

## takes two indices k and ell among [K]
## and returns feasible set corresponding to data N_k, N_ell, 0
f <- function(k, ell){
  countk <- counts[k]
  countell <- counts[ell]
  feasible_interval_kell <- sample_feasible_set(countk, countell)
  vertex1 <- rep(0, 3)
  vertex2 <- rep(0, 3)
  vertex3 <- rep(0, 3)
  vertex1[k] <- 1-feasible_interval_kell[1]
  vertex1[ell] <- feasible_interval_kell[1]
  vertex2[k] <- 1-feasible_interval_kell[2]
  vertex2[ell] <- feasible_interval_kell[2]
  vertex3[setdiff(1:3, c(k, ell))] <- 1
  vertex1_cart <- barycentric2cartesian(vertex1, v_cartesian)
  vertex2_cart <- barycentric2cartesian(vertex2, v_cartesian)
  vertex3_cart <- barycentric2cartesian(vertex3, v_cartesian)
  data.frame(x = c(vertex1_cart[1], vertex2_cart[1], vertex3_cart[1]), y = c(vertex1_cart[2], vertex2_cart[2], vertex3_cart[2]))
}

g <- gtriangle + geom_polygon(data = f(3,2), aes(x = x, y = y), alpha = 0.5)
g


k <- 2
ell <- 1
countk <- counts[k]
countell <- counts[ell]
feasible_interval_kell <- sample_feasible_set(countk, countell)
vertex1 <- rep(0, 3)
vertex2 <- rep(0, 3)
vertex3 <- rep(0, 3)
vertex1[k] <- 1-feasible_interval_kell[1]
vertex1[ell] <- feasible_interval_kell[1]
vertex2[k] <- 1-feasible_interval_kell[2]
vertex2[ell] <- feasible_interval_kell[2]
vertex3[setdiff(1:3, c(k, ell))] <- 1
vertex1_cart <- barycentric2cartesian(vertex1, v_cartesian)
vertex2_cart <- barycentric2cartesian(vertex2, v_cartesian)
vertex3_cart <- barycentric2cartesian(vertex3, v_cartesian)
polygon.df1 <- data.frame(x = c(vertex1_cart[1], vertex2_cart[1], vertex3_cart[1]), y = c(vertex1_cart[2], vertex2_cart[2], vertex3_cart[2]))


etas_ <- interval_to_etas(feasible_interval_kell)
etas <- diag(1, 3, 3)
etas[k,ell] <- etas_[1]
etas[ell,k] <- etas_[2]
etas[k, setdiff(1:3, c(k, ell))] <- Inf
etas[ell, setdiff(1:3, c(k, ell))] <- Inf
etas[setdiff(1:3, c(k, ell)), k] <- Inf
etas[setdiff(1:3, c(k, ell)), ell] <- Inf
cvx <- etas2cvxpolytope(etas)
#
cvx$vertices_barcoord
rbind(vertex1, vertex2, vertex3)

polygon.df2 <- data.frame(t(apply(cvx$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian))))
gtriangle + geom_polygon(data =polygon.df1, aes(x = x, y = y), alpha = 0.5, fill = "blue") +
  geom_polygon(data =polygon.df2, aes(x = X1, y = X2), alpha = 0.5, fill = "white", col = "black", linetype = 2) 



g <- function(counts, k, ell){
  countk <- counts[k]
  countell <- counts[ell]
  feasible_interval_kell <- sample_feasible_set(countk, countell)
  etas_ <- interval_to_etas(feasible_interval_kell)
  etas <- diag(1, 3, 3)
  etas[k,ell] <- etas_[1]
  etas[ell,k] <- etas_[2]
  etas[k, setdiff(1:3, c(k, ell))] <- Inf
  etas[ell, setdiff(1:3, c(k, ell))] <- Inf
  etas[setdiff(1:3, c(k, ell)), k] <- Inf
  etas[setdiff(1:3, c(k, ell)), ell] <- Inf
  return(etas)
}


##
## get three regions and intersect them
r1 <- g(counts, 1, 2)
cvx1 <- etas2cvxpolytope(r1)
polygon1.df <- data.frame(t(apply(cvx1$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian))))
r2 <- g(counts, 2, 3)
cvx2 <- etas2cvxpolytope(r2)
polygon2.df <- data.frame(t(apply(cvx2$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian))))
r3 <- g(counts, 1, 3)
cvx3 <- etas2cvxpolytope(r3)
polygon3.df <- data.frame(t(apply(cvx3$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian))))
##
cvx_intersection <- etas2cvxpolytope(pmin(r1, r2, r3))
cvx_intersection$vertices_barcoord

if (dim(cvx_intersection$vertices_barcoord)[1]==0){
  print("no intersection")
  cvxpolytope.df <- data.frame(X1 = 0, X2 = 0)
} else {
  cvxpolytope.cart <- t(apply(cvx_intersection$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  average_ <- colMeans(cvxpolytope.cart)
  o_ <- order(apply(sweep(cvxpolytope.cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvxpolytope.cart <- cvxpolytope.cart[o_,]
  cvxpolytope.df <- data.frame(cvxpolytope.cart)
  
}

gtriangle + geom_polygon(data = polygon1.df, aes(x = X1, y = X2), alpha = 0.5, fill = "blue") +
  geom_polygon(data = polygon2.df, aes(x = X1, y = X2), alpha = 0.5, fill = "red") +
  geom_polygon(data = polygon3.df, aes(x = X1, y = X2), alpha = 0.5, fill = "yellow") +
  geom_polygon(data = cvxpolytope.df, aes(x = X1, y = X2), alpha = 0.5, linetype = 2, col = "black")


# ### compute intersection of convex polytopes
# intersect_3polytopes <- function(cvx1, cvx2, cvx3){
#   test_constr <- cvx1$constr
#   test_constr$constr <- rbind(test_constr$constr, cvx2$constr$constr, cvx3$constr$constr)
#   test_constr$rhs <- c(test_constr$rhs, cvx2$constr$rhs, cvx3$constr$rhs)
#   test_constr$dir <- c(test_constr$dir, cvx2$constr$dir, cvx3$constr$dir)
#   ## make H representation
#   h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
#   h <- rcdd::redundant(h, representation = "H")$output
#   ## try to find V representation
#   v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
#   intersects <- (dim(v)[1] != 0)
#   return(list(v = v, h = h, test_constr = test_constr, intersects = intersects))
# }
# int_ <- intersect_3polytopes(cvx3, cvx1, cvx2)
# if (int_$intersects){
#   cvxpolytope <- int_$v[,-c(1,2)]
#   cvxpolytope.cart <- t(apply(t(apply(cvxpolytope, 1, function(v) c(v, 1-sum(v)))), 1, function(v) barycentric2cartesian(v, v_cartesian)))
#   average_ <- colMeans(cvxpolytope.cart)
#   o_ <- order(apply(sweep(cvxpolytope.cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
#   cvxpolytope.cart <- cvxpolytope.cart[o_,]
#   cvxpolytope.df <- data.frame(cvxpolytope.cart)
# } else {
#   cvxpolytope <- NULL
#   cvxpolytope.df <- data.frame(X1 = 0, X2 = 0)
# }


## reconstruct etas 


## new rejection sampler for K = 3
new_rejection_sampler <- function(counts){
  gcounts <- function(k,ell){
      countk <- counts[k]
      countell <- counts[ell]
      feasible_interval_kell <- sample_feasible_set(countk, countell)
      etas_ <- interval_to_etas(feasible_interval_kell)
      etas <- diag(1, 3, 3)
      etas[k,ell] <- etas_[1]
      etas[ell,k] <- etas_[2]
      etas[k, setdiff(1:3, c(k, ell))] <- Inf
      etas[ell, setdiff(1:3, c(k, ell))] <- Inf
      etas[setdiff(1:3, c(k, ell)), k] <- Inf
      etas[setdiff(1:3, c(k, ell)), ell] <- Inf
      return(etas)
  }
  K <- 3
  accept <- FALSE
  etas <- NULL
  nattempts <- 0
  while (!accept){
    r1 <- gcounts(1, 2)
    r2 <- gcounts(2, 3)
    r3 <- gcounts(3, 1)
    # etas <- pmin(r1, r2)
    etas <- pmin(r1, r2, r3)
    nattempts <- nattempts + 1
    accept <- dempsterpolytope:::check_cst_graph(etas)
    ## alternative: get convex polytope
    # cvx_intersection <- etas2cvxpolytope(etas)
    ## and check whether dim(cvx_intersection$vertices_barcoord)[1] == 0
  }
  return(list(etas = etas, nattempts = nattempts))
}


dempsterpolytope:::rejectionsampler(counts)
new_rejection_sampler(counts)

nsamples_rs <- 5e2
samples_rs <- foreach(irep = 1:nsamples_rs) %dorng% {
  dempsterpolytope::rejectionsampler(counts)
}

samples_rs2 <- foreach(irep = 1:nsamples_rs) %dorng% {
  new_rejection_sampler(counts)
}

library(plyr)
res <- laply(samples_rs, function(l) as.matrix(l$etas))
res2 <- laply(samples_rs2, function(l) as.matrix(l$etas))

param <- 3
xgrid <- seq(from = 0, to = 1, length.out = 100)
cdf_rs <- etas_to_lower_upper_cdf_dopar(res, param, xgrid)
lowercdf_rs <- colMeans(cdf_rs$iscontained)
uppercdf_rs <- colMeans(cdf_rs$intersects)
cdf_rs2 <- etas_to_lower_upper_cdf_dopar(res2, param, xgrid)
lowercdf_rs2 <- colMeans(cdf_rs2$iscontained)
uppercdf_rs2 <- colMeans(cdf_rs2$intersects)
#
plot(xgrid, lowercdf_rs, type = "l")
lines(xgrid, uppercdf_rs)
lines(xgrid, lowercdf_rs2, col = "blue", lty = 2)
lines(xgrid, uppercdf_rs2, col = "blue", lty = 2)


