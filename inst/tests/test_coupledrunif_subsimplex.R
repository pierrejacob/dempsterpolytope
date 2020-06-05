## Test of various couplings of uniform sampling on Delta_k(theta), 
## the subsimplex where vertex k is replaced by theta
## 
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

n <- 20000
K <- 4
theta_star1 <- rexp(K)
theta_star1 <- theta_star1 / sum(theta_star1)
theta_star2 <- rexp(K)
theta_star2 <- theta_star2 / sum(theta_star2)
k <- 1

## independent sampling
pts1 <- dempsterpolytope:::runif_piktheta_cpp(n, k, theta_star1)$pts
mean(log(pts1[,1]) * pts1[,K]^2)
pts2 <- dempsterpolytope:::runif_piktheta_cpp(n, k, theta_star2)$pts
mean(log(pts2[,1]) * pts2[,K]^2)

## common random numbers
pts_ <- dempsterpolytope:::crng_runif_piktheta_cpp(n, k, theta_star1, theta_star2)
mean(log(pts_$pts1[,1]) * pts_$pts1[,K]^2)
mean(log(pts_$pts2[,1]) * pts_$pts2[,K]^2)

rmaxcoupling <- function(k, theta_star1, theta_star2, multiplier = 0.9){
  x <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star1)
  pdf1_x <- 1/theta_star1[k]
  pdf2_x <- dempsterpolytope:::dunif_piktheta_cpp(x, k, theta_star2)
  if (runif(1) < (multiplier * pdf2_x / pdf1_x)){
    return(list(pts = cbind(x,x), equal = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star2)
      pdf2_y <- 1/theta_star2[k]
      pdf1_y <- dempsterpolytope:::dunif_piktheta_cpp(y, k, theta_star1)
      reject <- (runif(1) < min(multiplier, pdf1_y/pdf2_y))
    }
    return(list(pts = cbind(x,y), equal = FALSE))
  }
}

print(theta_star1)
pts_maxcoupled <- foreach(irep = 1:n) %dorng% {
  rmaxcoupling(k, theta_star1, theta_star2)
}
pts_maxcoupled1 <- t(sapply(pts_maxcoupled, function(x) x$pts[,1]))
pts_maxcoupled2 <- t(sapply(pts_maxcoupled, function(x) x$pts[,2]))

mean(log(pts_maxcoupled1[,1]) * pts_maxcoupled1[,K]^2)
mean(log(pts_maxcoupled2[,1]) * pts_maxcoupled2[,K]^2)

mean(sapply(pts_maxcoupled, function(x) x$equal))
# 
# print(theta_star1)
pts_maxcoupled_cpp <- dempsterpolytope:::maxcoupling_runif_piktheta_cpp(n, k, theta_star1, theta_star2)
mean(log(pts_maxcoupled_cpp$pts1[,1]) * pts_maxcoupled_cpp$pts1[,K]^2)
mean(log(pts_maxcoupled_cpp$pts2[,1]) * pts_maxcoupled_cpp$pts2[,K]^2)


pts_maxcoupled_cpp$ncoupled / n
# pts_maxcoupled_cpp$minratios1
# pts_maxcoupled_cpp$minratios2

##


