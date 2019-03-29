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

## plot points of the form 1/2 + phi/4, 1/2 - phi/2, phi/4

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

## generate counts
phi_star <- 0.5
theta_star <- A * phi_star  + b
theta_star_cart <- barycentric2cartesian(theta_star, v_cartesian)
n <- 6
X <- sample(x = 1:K, size = n, replace = TRUE, prob = theta_star)
freqX <- tabulate(X)
freqX
##
pts <- initialize_pts(freqX, theta_star)
etas <- do.call(rbind, pts$minratios)

pts_cart <- lapply(pts$pts, function(l) t(apply(l, 1, function(v) barycentric2cartesian(v, v_cartesian))))
##
g <- ggplot_triangle(v_cartesian, pts_cartesian = pts_cart) 
g <- add_segment_linkage(g, 0, 1, A, b, v_cartesian, colour = "black")
g <- g + geom_point(aes(x = theta_star_cart[1], y = theta_star_cart[2]), size = 3, colour = 'orange')
g <- g + geom_polygon(data = data.frame(x = c(theta_star_cart[1], v2[1], v3[1]), y = c(theta_star_cart[2], v2[2], v3[2])),
                      aes(fill = "1"), colour = "black", alpha = 0.15)
g <- g + geom_polygon(data = data.frame(x = c(theta_star_cart[1], v1[1], v3[1]), y = c(theta_star_cart[2], v1[2], v3[2])),
                      aes(fill = "2"), colour = "black", alpha = 0.15)
g <- g + geom_polygon(data = data.frame(x = c(theta_star_cart[1], v1[1], v2[1]), y = c(theta_star_cart[2], v1[2], v2[2])),
                      aes(fill = "3"), colour = "black", alpha = 0.15)
g <- g + theme(legend.position = "none")
g
## 

## eta_{j, i} is etas[j,i]
## get lower and upper bound on phi from the ratios
get_lower_upper <- function(etas, A, b){
  upper <- 1
  lower <- 0
  K <- dim(etas)[1]
  for (k1 in 1:K){
    for (k2 in setdiff(1:K, k1)){
      c <- etas[k1,k2]
      denom <- A[k2] - c * A[k1]
      if (denom >= 0){
        upper <- min(upper, (c * b[k1] - b[k2]) / (denom))
      } else {
        lower <- max(lower, (c * b[k1] - b[k2]) / (denom))
      }
    }
  }
  return(c(lower, upper))
}
lu <- get_lower_upper(etas, A, b)
lower <- lu[1]; upper <- lu[2]
g <- ggplot_triangle(v_cartesian, pts_cartesian = pts_cart, etas = etas, removelines = FALSE, addpolytope = T) 
g <- add_segment_linkage(g, lower, upper, A, b, v_cartesian, colour = "yellow")
g
## 
## Now refresh categories
## get lower and upper bound on 
get_lower_upper_updatek <- function(etas, A, b, k){
  upper <- 1
  lower <- 0
  K <- dim(etas)[1]
  for (k1 in setdiff(1:K, k)){
    for (k2 in setdiff(1:K, k1)){
      c <- etas[k1,k2]
      denom <- A[k2] - c * A[k1]
      if (denom >= 0){
        upper <- min(upper, (c * b[k1] - b[k2]) / (denom))
      } else {
        lower <- max(lower, (c * b[k1] - b[k2]) / (denom))
      }
    }
  }
  return(c(lower, upper))
}
k <- 1
lu_k <- get_lower_upper_updatek(etas, A, b, k)
# now compute the corresponding thetas, and get the corresponding k-th component
(A * lu_k[1] + b)[k]
(A * lu_k[2] + b)[k]
# and get the phi corresponding to maximal theta_k
if ((A * lu_k[1] + b)[k] < (A * lu_k[2] + b)[k]){
  # then max theta_k correspond to phi = upper bound
  theta_star <- A * lu_k[2] + b
} else {
  theta_star <- A * lu_k[1] + b
}
##

gibbs_sampler_linkage <- function(niterations, freqX, phi_0, A, b){
  K_ <- length(freqX)
  if (missing(phi_0)){
    phi_0 <- 0.5
  }
  theta_0 <- A * phi_0 + b
  categories <- 1:K_
  # store points in barycentric coordinates
  Achain <- list()
  for (k in 1:K_){
    if (freqX[k] > 0){
      Achain[[k]] <- array(0, dim = c(niterations, freqX[k], K_))
    } else {
      Achain[[k]] <- array(0, dim = c(niterations, 1, K_))
    }
  }
  # store constraints in barycentric coordinates
  etas_chain <- array(0, dim = c(niterations, K_, K_))
  lu_chain <- matrix(0, nrow = niterations, ncol = 2)
  ## initialization
  init_tmp <- initialize_pts(freqX, theta_0)
  pts <- init_tmp$pts
  # store points
  for (k in 1:K_){
    if (freqX[k] > 0){
      Achain[[k]][1,,] <- pts[[k]]
    } else {
      Achain[[k]][1,1,] <- rep(1/(K_-1), K_)
      Achain[[k]][1,1,k] <- 0
    }
  }
  etas <- do.call(rbind, init_tmp$minratios)
  # store constraints
  etas_chain[1,,] <- etas
  lu_chain[1,] <- get_lower_upper(etas, A, b)
  # loop over Gibbs sampler iterations
  for (iter_gibbs in 2:niterations){
    # loop over categories
    for (k in categories){
      if (freqX[k] > 0){
        lu_k <- get_lower_upper_updatek(etas, A, b, k)
        if ((A * lu_k[1] + b)[k] < (A * lu_k[2] + b)[k]){
          # then max theta_k correspond to phi = upper bound
          theta_star <- A * lu_k[2] + b
        } else {
          # then max theta_k correspond to phi = lower bound
          theta_star <- A * lu_k[1] + b
        }
        ##
        pts_k <- montecarlodsm:::runif_piktheta_cpp(freqX[k], k, theta_star)
        pts[[k]] <- pts_k$pts
        etas[k,] <- pts_k$minratios
      }
    }
    # store points and constraints
    for (k in categories){
      if (freqX[k] > 0){
        Achain[[k]][iter_gibbs,,] <- pts[[k]]
      } else {
        Achain[[k]][iter_gibbs,1,] <- rep(1/(K_-1), K_)
        Achain[[k]][iter_gibbs,1,k] <- 0
      }
    }
    etas_chain[iter_gibbs,,] <- etas
    lu_chain[iter_gibbs,] <- get_lower_upper(etas, A, b)
  }
  # return points post-burnin
  return(list(etas_chain = etas_chain, Achain = Achain, lu_chain = lu_chain))
}

gibbs_res_ <- gibbs_sampler_linkage(100, freqX, A = A, b = b)
# gibbs_res_$etas_chain[100,,]
# gibbs_res_$lu_chain[100,]
# all(apply(gibbs_res_$lu_chain, 1, function(v) v[1] < v[2]))
# all(apply(gibbs_res_$lu_chain, 1, function(v) 0 < v[1]))
# all(apply(gibbs_res_$lu_chain, 1, function(v) v[2] < 1))

g <- ggplot_triangle(v_cartesian, etas = gibbs_res_$etas_chain[100,,], addpolytope = TRUE)
g <- add_segment_linkage(g, gibbs_res_$lu_chain[100,1], gibbs_res_$lu_chain[100,2], A, b, v_cartesian, colour = "yellow")
g

## now retrieve segment of phi from polytope on theta
etas <- gibbs_res_$etas_chain[100,,]
# construct polytope
cvx <- etas2cvxpolytope(etas)
# now construct linear constraints defining the phi-theta relationship
# theta = A * phi + b so <= and >=
# A_1 is positive so theta_1 is in (b_1, A_1 + b_1)
# and then theta_j for j > 1 is deterministic given theta_1:
# writing phi as (theta_1 - b_1)/ A_1 
# theta_j = A_j * (theta_1 - b_1)/ A_1 + b_j 
# i.e. 
# -A_j/A_i * theta_1 + theta_j = b_j - b_1 * A_j/A_1
constr <- matrix(0, nrow = 2, ncol = K)
rhs <- c(-b[1], A[1] + b[1])
constr[1,1] <- -1
constr[2,1] <- 1
dir <- c("<=", "<=")
for (j in 2:K){
  row_constr_ <- rep(0, K)
  row_constr_[1] <- -A[j]/A[1]
  row_constr_[j] <- 1
  rhs <- c(rhs, b[j] - b[1] * A[j] / A[1], -(b[j] - b[1] * A[j] / A[1]))
  constr <- rbind(constr, row_constr_, -row_constr_, deparse.level = 0)
  dir <- c(dir, "<=", "<=")
}

# check that it works
hitandrun::findVertices(list(constr = constr, dir = dir, rhs = rhs))
A * 0 + b
A * 1 + b

## can we plot the lower and upper CDF of phi now?
K <- 4
freqX <- c(25, 3, 4, 7)
A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)
gibbs_res_2 <- gibbs_sampler_linkage(10000, freqX, A = A, b = b)
# all(apply(gibbs_res_2$lu_chain, 1, function(v) v[1] < v[2]))
# all(apply(gibbs_res_2$lu_chain, 1, function(v) 0 < v[1]))
# all(apply(gibbs_res_2$lu_chain, 1, function(v) v[2] < 1))
lu_chain <- gibbs_res_2$lu_chain[100:(dim(gibbs_res_2$lu_chain)[1]),]
matplot(lu_chain, type = "l")
acf(lu_chain[,1])
acf(lu_chain[,2])
# for a given phi_0, and assertion {phi < phi_0}
# upper probability is proportion of time random interval [phi_lower, phi_upper] intersects with [0, phi_0]
# lower probability is proportion of time random interval [phi_lower, phi_upper] is contained in [0, phi_0]
cdf_lowerupper <- function(phi_0, lu_chain){
  # intersect if phi_lower < phi_0
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < phi_0))
  # contains if phi_upper < phi_0 
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < phi_0))
  return(c(cdf_lower, cdf_upper))
}

phi_grid_length <- 200
phi_grid <- seq(from = 0, to = 1, length.out = phi_grid_length)
cdf_lower_values <- rep(0, phi_grid_length)
cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- cdf_lowerupper(phi_grid[iphi], lu_chain)
  cdf_lower_values[iphi] <- cdf_res[1]
  cdf_upper_values[iphi] <- cdf_res[2]
}
plot(phi_grid, cdf_lower_values, type = "l")
lines(phi_grid, cdf_upper_values)
# there you go, we can compare with Figure 6.2 in Lawrence et al.
