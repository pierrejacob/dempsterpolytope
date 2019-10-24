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
n <- 20
## number of categories
K <- 3
categories <- 1:K
## data 
counts <- c(7,5,8)
##
## encompassing triangle has three vertices
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")

###
niterations <- 200
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts, theta_0 = c(0.8,0.1,0.1))
(proc.time() - pct)[3]
##


## now get all polytopes of feasible parameters at all iterations
## and overlay them in plot
df.polytope <- data.frame()
for (iteration in 1:niterations){
  etas_iteration <- samples_gibbs$etas[iteration,,]
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

### Show an instance of a feasible plot with all six linear constraints
iteration <- 170
etas <- samples_gibbs$etas[iteration,,]
pts_barcoord <- lapply(samples_gibbs$Us, function(l) l[iteration,,])
pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))
g <- ggplot_triangle(v_cartesian, pts_cart, etas, addpolytope = T, cols = cols)
g

