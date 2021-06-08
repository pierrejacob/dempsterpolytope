## This scripts runs the Gibbs sampler 
## on counts in 3 categories
## and shows resulting polytopes in plots

rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
set.seed(4)

## number of observations
n <- 20
## number of categories
K <- 3
categories <- 1:K
## data 
counts <- c(7,5,8)
###
niterations <- 200
pct <- proc.time()
gibbs_results <- gibbs_sampler(niterations = niterations, counts = counts)
(proc.time() - pct)[3]
##

## now get vertices for all polytopes 
## and overlay them in a plot
df.polytope <- data.frame()
for (iteration in 1:niterations){
  etas_iteration <- gibbs_results$etas[iteration,,]
  vertices <- etas_vertices(etas_iteration)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(vertices, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                   iteration = iteration))
}
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=df.polytope %>% filter(iteration >= 100), aes(x = x, y = y, group = iteration), alpha = .2,
               colour = 'black', size = 0.25)
g


