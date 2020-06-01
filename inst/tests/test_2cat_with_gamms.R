library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())
##
## K = 2
v_cartesian <- list(c(0, 0), c(1,0))
cols <- c("red", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]
## Data
counts <- c(2,5)
N <- sum(counts)
## Gibbs sampler
burnin <- 1000
niterations <- 2000
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)

iter <- sample.int(niterations, 1)
cvxp <- etas2cvxpolytope(samples_gibbs$etas[iter,,])
cvxp

cvxp$vertices_cart <- t(apply(cvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
plot(x = c(0,1), y = c(0,0), type = "l")
abline(v = cvxp$vertices_cart[,1], lty = 2)

## otherwise from the etas

## theta_1/theta_2 <= eta_2->1 
## theta_2/theta_1 <= eta_1->2 
## implies 1 / (1 + eta_{1->2}) <= theta_1 <= eta_{2->1} / (1 + eta_{2->1})
c(1 - 1/(1+samples_gibbs$etas[iter,1,2]), 1 - samples_gibbs$etas[iter,2,1]/(1+samples_gibbs$etas[iter,2,1]))
## which exactly matches
cvxp$vertices_cart[,1]

param <- 1
nsubiterations <- 1e4
subiterations <- floor(seq(from = burnin+1, to = niterations, length.out = nsubiterations))
xgrid <- seq(from = 0, to = 1, length.out = 100)
cdfs_ <- etas_to_lower_upper_cdf_dopar(samples_gibbs$etas[subiterations,,], param, xgrid)
lowercdf <- colMeans(cdfs_$iscontained)
uppercdf <- colMeans(cdfs_$intersects)

# lower / upper
plot(xgrid, lowercdf, type = "l", xlab = "theta1", ylab = "CDF")
lines(xgrid, uppercdf, lty = 2)
##

## exact sampling using sorted uniforms
sample_feasible_set <- function(counts){
  unsorted_unifs <- runif(N)
  sorted_unifs <- sort(unsorted_unifs)
  interval_ <- sorted_unifs[c(counts[1], counts[1]+1)]
  return(interval_)
}

sample_feasible_set2 <- function(counts){
  gamma1 <- rgamma(n = 1, shape = counts[1], rate = 1)
  gamma2 <- rgamma(n = 1, shape = counts[2], rate = 1)
  wi <- rexp(1, 1)
  interval_ <- c(gamma1/(gamma1+wi+gamma2), (gamma1+wi)/(gamma1+wi+gamma2))
  return(interval_)
}
## feasible sets
feasible_sets_ <- foreach(irep = 1:1e4, .combine = rbind) %dorng% sample_feasible_set2(counts)
# feasible_sets_ %>% head
## lowercdf: proportion of feasible sets that are contained in [0,x], i.e. upper end < x
## uppercdf: proportion of feasible sets that intersect in [0,x] i.e. lower end < x
alt_lowercdf <- rep(0, length(xgrid))
alt_uppercdf <- rep(0, length(xgrid))
for (igrid in 1:length(xgrid)){
  x <- xgrid[igrid]
  alt_lowercdf[igrid] <- mean(apply(feasible_sets_, 1, function(v) v[2] < x))
  alt_uppercdf[igrid] <- mean(apply(feasible_sets_, 1, function(v) v[1] < x))
}

plot(xgrid, lowercdf, type = "l", xlab = "theta1", ylab = "CDF")
lines(xgrid, uppercdf, lty = 2)
lines(xgrid, alt_lowercdf, lty = 3, col = "blue")
lines(xgrid, alt_uppercdf, lty = 3, col = "blue")

####
dim(samples_gibbs$Us[[1]])
dim(samples_gibbs$Us[[2]])

##
iter <- sample.int(niterations, 1)
u1 <- samples_gibbs$Us[[1]][iter,,]
u2 <- samples_gibbs$Us[[2]][iter,,]
feasible_interval <- c(1 - 1/(1+samples_gibbs$etas[iter,1,2]), 1 - samples_gibbs$etas[iter,2,1]/(1+samples_gibbs$etas[iter,2,1]))

## plot points in segment 
plot(x = c(0,1), y = c(0,0), type = "l", xlab = "", ylab = "", ylim = c(-.1, 1.1))
## 2 points in category 1, close to vertex 2
points(x = u1[,2], y = rep(0, nrow(u1)), pch = 2)
## 5 points in category 2, close to vertex 1
points(x = u2[,2], y = rep(0, nrow(u2)), pch = 1)
abline(v = feasible_interval, lty = 3)

## add empty category
## extend us to add third component
extend_to_K3 <- function(v){
  sum_ <- rgamma(1, 2, 1)
  y <- c(v * sum_, rexp(1, 1))
  return(y / sum(y))
}
u1_extended <- t(apply(u1, 1, extend_to_K3))
u2_extended <- t(apply(u2, 1, extend_to_K3))


# 
v_cartesian <- list(c(0,0), c(1,0), c(1/2, sin(pi/3)))
# # cols <- c("red", "yellow", "blue")
# # contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))

u1_extended_cart <- t(apply(u1_extended, 1, function(v) barycentric2cartesian(v, v_cartesian)))
u2_extended_cart <- t(apply(u2_extended, 1, function(v) barycentric2cartesian(v, v_cartesian)))

g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_point(data = data.frame(u1_extended_cart), aes(x = X1, y = X2), pch = 2)
g <- g + geom_point(data = data.frame(x = u1[,2], y = rep(0, nrow(u1))), aes(x = x, y = y), pch = 2, col = "red")
g <- g + geom_segment(data = data.frame(x = u1[,2], y = rep(0, nrow(u1)), xend = u1_extended_cart[,1], yend = u1_extended_cart[,2]), aes(x = x, y = y, xend = xend, yend = yend), col = "red")
g <- g + geom_point(data = data.frame(u2_extended_cart), aes(x = X1, y = X2), pch = 1)
g <- g + geom_point(data = data.frame(x = u2[,2], y = rep(0, nrow(u2))), aes(x = x, y = y), pch = 1, col = "red")
g <- g + geom_segment(data = data.frame(x = u2[,2], y = rep(0, nrow(u2)), xend = u2_extended_cart[,1], yend = u2_extended_cart[,2]), aes(x = x, y = y, xend = xend, yend = yend), col = "red")
## feasible set
feasible_set_K3 <- data.frame(x = c(sort(feasible_interval)[1], sort(feasible_interval)[2], v3[1]), y = c(0,0,v3[2]))
g <- g + geom_polygon(data = feasible_set_K3, aes(x = x, y = y), alpha = 0.5)
g

## retrieve feasible set using etas
etas <- diag(1, 3, 3)
etas[1,2] <- min(u1_extended[,2]/u1_extended[,1])
etas[2,1] <- min(u2_extended[,1]/u2_extended[,2])
# min(u1_extended[,2]/u1_extended[,1]) == min(u1[,2]/u1[,1])
etas[3,1] <- etas[3,2] <- Inf
etas[1,3] <- etas[2,3] <- Inf
etas
feasible_set_cvxp <- etas2cvxpolytope(etas)
vertices_cart <- t(apply(feasible_set_cvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
# order vertices according to angles
average_ <- colMeans(vertices_cart)
o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
vertices_cart <- vertices_cart[o_,]
df.polytope <- data.frame(x = vertices_cart[,1], y= vertices_cart[,2])
g + geom_polygon(data = df.polytope, aes(x = x, y = y), alpha = 0.5, fill = "blue")


# g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)


# 
# ## data with K = 3 categories
# counts <- c(2,5)
# ## run Gibbs sampler
# niterations <- 500
# samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
# ##
# # dim(samples_gibbs$etas)
# # samples_gibbs$etas[3,,]
# 
# ## overlay feasible regions in simplex
# 
# iteration <- 10
# etas <- samples_gibbs$etas[iteration,,]
# pts_barcoord <- lapply(samples_gibbs$Us, function(l) l[iteration,,])
# pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))
# ##
# pts_barcoord
# pts_cart
# 
# g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
# g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
# g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
# g
