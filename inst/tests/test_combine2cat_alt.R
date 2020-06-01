## This implements Art's idea
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
counts <- c(3,1,3)
N <- sum(counts)
new_rejection_sampler <- function(counts){
  ## draw Gamma for each vertex
  gamma_vertices <- sapply(counts, function(shape) rgamma(n = 1, shape = shape, rate = 1))
  names(gamma_vertices) <- 1:K
  exp_edges <- c()
  etas <- diag(1, 3, 3)
  for (k in 1:(K-1)){
    for (ell in setdiff(k:K, k)){
      gamma_k <- gamma_vertices[paste0(k)]
      gamma_ell <- gamma_vertices[paste0(ell)]
      exp_ <- rexp(n = 1, rate = 1)
      interval <- c(gamma_k/(gamma_k+exp_+gamma_ell), (gamma_k+exp_)/(gamma_k+exp_+gamma_ell))
      etaktoell <- interval[2] / (1 - interval[2])
      etaelltok <- (1 - interval[1]) / interval[1]
      etas[k,ell] <- etaktoell
      etas[ell,k] <- etaelltok
    }
  }
  return(list(etas = etas))
}
cvx <- etas2cvxpolytope(new_rejection_sampler(counts)$etas)
#
cvx$vertices_barcoord
cvxpolytope.cart <- t(apply(cvx$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
average_ <- colMeans(cvxpolytope.cart)
o_ <- order(apply(sweep(cvxpolytope.cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
cvxpolytope.cart <- cvxpolytope.cart[o_,]
cvxpolytope.df <- data.frame(cvxpolytope.cart)

g <- gtriangle + geom_polygon(data = cvxpolytope.df, aes(x = X1, y = X2), alpha = 0.5)
g

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

param <- 2
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


