library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(2)
rm(list = ls())

## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
##
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
#
A <- c(1/4, -1/2, 1/4)
b <- c(1/2,  1/2,  0)
# theta = A phi + b, phi in (0,1)

# Now with ggplot2
g <- ggplot_triangle(v_cartesian)
add_segment_linkage <- function(g, lower, upper, A, b, v_cartesian, colour = "yellow"){
  theta_min <- lower * A + b
  theta_max <- upper * A + b
  theta_min_cart <- barycentric2cartesian(theta_min, v_cartesian)
  theta_max_cart <- barycentric2cartesian(theta_max, v_cartesian)
  return(g + geom_segment(data = data.frame(x = theta_min_cart[1], xend = theta_max_cart[1], y = theta_min_cart[2], yend = theta_max_cart[2]), 
                          aes(x = x, y = y, xend = xend, yend = yend), colour = colour))
} 

g <- add_segment_linkage(g, 0, 1, A, b, v_cartesian, colour = "black")
print(g)


### run standard Gibbs sampler (without linkage constrained)
phi_star <- 0.5
theta_star <- A * phi_star  + b
n <- 1000
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_star)
freqX <- tabulate(X)
freqX

niterations <- 2e4
warmup <- 1e3
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
##

nsubiterations <- 5000
subiterations <- floor(seq(from = warmup, to = niterations, length.out = nsubiterations))
## and overlay them in plot
df.polytope <- data.frame()
for (index in 1:nsubiterations){
  etas <- samples_gibbs$etas_chain[subiterations[index],,]
  etascvxp <- etas2cvxpolytope(etas)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                               iteration = index))
}
g <- ggplot_triangle(v_cartesian) +
  geom_polygon(data=df.polytope, aes(x = x, y = y, group = iteration), alpha = 1)
g <- add_segment_linkage(g, 0, 1, A, b, v_cartesian, colour = "black")
print(g)

## add theta DGP
theta_star_cart <- barycentric2cartesian(theta_star, v_cartesian)
theta_mle_cart <- barycentric2cartesian(freqX/n, v_cartesian)
g + geom_point(aes(x = theta_star_cart[1], y = theta_star_cart[2])) + geom_point(aes(x = theta_mle_cart[1], y = theta_mle_cart[2]), col = "red")

# ## sub-set etas which intersect with "linkage" constraints
# etas <- samples_gibbs$etas_chain[subiterations,,]
# 
# eta <- etas[10,,]
# cvxp <- etas2cvxpolytope(eta)
# ## visualize that particular eta
# vertices_cart <- t(apply(cvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
# average_ <- colMeans(vertices_cart)
# o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
# vertices_cart <- vertices_cart[o_,]
# df.polytope <- data.frame(x = vertices_cart[,1], y= vertices_cart[,2])
# 
# g + geom_polygon(data = df.polytope, fill = "yellow")


# now construct linear constraints defining the phi-theta relationship
# theta = A * phi + b so <= and >=
# A_1 is positive so theta_1 is in (b_1, A_1 + b_1)
# and then theta_j for j > 1 is deterministic given theta_1:
# writing phi as (theta_1 - b_1)/ A_1 
# theta_j = A_j * (theta_1 - b_1)/ A_1 + b_j 
# i.e. 
# -A_j/A_i * theta_1 + theta_j = b_j - b_1 * A_j/A_1

constr <- matrix(0, nrow = 2, ncol = K-1)
rhs <- c(-b[1], A[1] + b[1])
constr[1,1] <- -1
constr[2,1] <- 1
dir <- c("<=", "<=")
for (j in 2:(K-1)){
  row_constr_ <- rep(0, K-1)
  row_constr_[1] <- -A[j]/A[1]
  row_constr_[j] <- 1
  rhs <- c(rhs, b[j] - b[1] * A[j] / A[1], -(b[j] - b[1] * A[j] / A[1]))
  constr <- rbind(constr, row_constr_, -row_constr_, deparse.level = 0)
  dir <- c(dir, "<=", "<=")
}
segment_constr <- list(constr = constr, dir = dir, rhs = rhs)
## check that it works
# hitandrun::findVertices(segment_constr)
# A * 0 + b
# A * 1 + b

intersects_with_segment <- function(eta, segment_constr){
  cvxp <- etas2cvxpolytope(eta)
  test_constr <- segment_constr
  test_constr$constr <- rbind(test_constr$constr, cvxp$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp$constr$dir)
  ## make H representation (H for ?)
  h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
  ## try to find V representation (V for Vendetta or Vertex?)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  intersects <- (dim(v)[1] != 0)
  return(intersects)  
}

etas <- samples_gibbs$etas_chain[subiterations,,]
intersects_segment_ <- sapply(1:dim(etas)[1], function(i) intersects_with_segment(etas[i,,], segment_constr))
sum(intersects_segment_)

df.polytope <- data.frame()
for (index in 1:nsubiterations){
  etas <- samples_gibbs$etas_chain[subiterations[index],,]
  if (intersects_with_segment(etas, segment_constr)){
    etascvxp <- etas2cvxpolytope(etas)
    ## convert coordinates to cartesian
    vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
    # order vertices according to angles
    average_ <- colMeans(vertices_cart)
    o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
    vertices_cart <- vertices_cart[o_,]
    df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                                 iteration = index))
  }
}
g <- ggplot_triangle(v_cartesian) +
  geom_polygon(data=df.polytope, aes(x = x, y = y, group = iteration), alpha = 1)
g <- add_segment_linkage(g, 0, 1, A, b, v_cartesian, colour = "black")
g + geom_point(aes(x = theta_star_cart[1], y = theta_star_cart[2])) + geom_point(aes(x = theta_mle_cart[1], y = theta_mle_cart[2]), col = "red")

