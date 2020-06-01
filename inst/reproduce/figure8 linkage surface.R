rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set.seed(1)
## 

x <- seq(-1,1,0.01) 
y <- seq(-1,1,0.01)
f <- function(x,y){ z <- -x - y + 1 }
z <- outer(x,y,f)
z <- ifelse(z<0,NA,z)
# persp_res <- persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "white", border = NA, box = FALSE, axes = 0)

v1 <- c(sqrt(8/9),0,-1/3)
v2 <- c(-sqrt(2/9),sqrt(2/3),-1/3)
v3 <- c(-sqrt(2/9),-sqrt(2/3),-1/3)
v4 <- c(0,0,1)
z <- matrix(rep(-1/3, 4), ncol = 2)
z[1,] <- 1
# persp_res <- persp(x = c(-sqrt(2/9), sqrt(8/9)), y = c(-sqrt(2/3), sqrt(2/3)), z = z, 
#                    theta = 200, phi = 1, expand = 1, col = "white", border = NA, box = FALSE, axes = 0)
persp_res <- persp(x = c(-sqrt(2/9), sqrt(8/9)), y = c(-sqrt(2/3), sqrt(2/3)), z = z, 
                   r = 2, d = 2,
                   xlim = c(-.5,.6), zlim = c(-.45,0.95),
                   theta = 200, phi = 1, expand = 1, col = "white", border = NA, box = FALSE, axes = 0)


add_segment <- function(v, vprime, persp_res){
  seg_ <- cbind(v, vprime)
  lines(trans3d(seg_[1,], seg_[2,], seg_[3,], pmat = persp_res), col = rgb(0,0,0))
}


add_segment(v1, v2, persp_res)
add_segment(v1, v3, persp_res)
add_segment(v2, v3, persp_res)
add_segment(v1, v4, persp_res)
add_segment(v2, v4, persp_res)
add_segment(v3, v4, persp_res)

A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)
# u_grid <- rev(seq(from = 0, to = 1, by = 0.02))
# for (phi in u_grid){
#   point_ <- (A * phi + b)[1] * v1 + (A * phi + b)[2] * v2 + (A * phi + b)[3] * v3 + (A * phi + b)[4] * v4 
#   points(trans3d(point_[1], point_[2], point_[3], pmat = persp_res), col = rgb(0,0,0))
# }
phi <- 0
point_0 <- (A * phi + b)[1] * v1 + (A * phi + b)[2] * v2 + (A * phi + b)[3] * v3 + (A * phi + b)[4] * v4
phi <- 1
point_1 <- (A * phi + b)[1] * v1 + (A * phi + b)[2] * v2 + (A * phi + b)[3] * v3 + (A * phi + b)[4] * v4 
seg_ <- cbind(point_0, point_1)
lines(trans3d(seg_[1,], seg_[2,], seg_[3,], pmat = persp_res), col = rgb(0,0,0), lty = 2)
points(trans3d(seg_[1,], seg_[2,], seg_[3,], pmat = persp_res), col = rgb(0,0,0))


text(trans3d(v1[1]+0.1,v1[2],v1[3],persp_res), "1", cex = 2)
text(trans3d(v2[1]-0.1,v2[2],v2[3],persp_res), "2", cex = 2)
text(trans3d(v3[1]-0.1,v3[2],v3[3],persp_res), "3", cex = 2)
text(trans3d(v4[1],v4[2],v4[3]+0.1,persp_res), "4", cex = 2)


###
# dev.copy(png,'myplot.png')
# dev.off()
dev.copy(pdf,"linkagesurface.pdf")
dev.off()

# text(trans3d(v1[1],v1[2],v1[3],persp_res), "1", cex = 1)
# text(trans3d(v2[1],v2[2],v2[3],persp_res), "2", cex = 1)
# text(trans3d(v3[1],v3[2],v3[3],persp_res), "3", cex = 1)
# text(trans3d(v4[1],v4[2],v4[3],persp_res), "4", cex = 1)

###

# define A and b such that theta = A phi + b where phi is in the interval (0,1)
A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)

# data 
K <- 4
counts = c(25, 3, 4, 7)

## let's choose a burn-in of 50, based on the above plot
burnin <- 50
nchains <- 25
niterations <- burnin+1
library(abind)
etas <- foreach(irep = 1:nchains) %dorng% {
  init <- rexp(K)
  init <- init/sum(init)
  samples_gibbs <- gibbs_sampler(niterations, counts, theta_0 = init)
  samples_gibbs$etas[(burnin+1):niterations,,]
}
dim(etas[[1]])
# barypoints <- etas2cvxpolytope(etas[[1]])$vertices_barcoord
## 
# cartpoints <- apply(barypoints, 1, function(x) x[1] * v1 + x[2] * v2 + x[3] * v3 + x[4] * v4)

# points(trans3d(cartpoints[1,], cartpoints[2,], cartpoints[3,], pmat = persp_res), col = rgb(0,0,0))


