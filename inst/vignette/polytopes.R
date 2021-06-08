## This scripts illustrates some manipulation of polytopes
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
set.seed(6)
## number of categories
K <- 3
categories <- 1:K
##
## data 
counts <- c(3,1,2)
## run Gibbs sampler
niterations <- 70
pct <- proc.time()
gibbs_results <- gibbs_sampler(niterations = niterations, counts = counts)
(proc.time() - pct)[3]
##

## get a KxK matrix of etas
etas <- gibbs_results$etas[niterations,,]
etas
## obtain the vertices of the convex polytope by enumeration
## where theta is in the simplex and 
## theta_ell / theta_k <= eta[k,l] for all k,l in [K]
vertices <- etas_vertices(etas)
## each row represents a vertex, and sums to one

## to plot vertices, we convert the barycentric coordinates to cartesian coordinates
# ?barycentric2cartesian
vertices_cartesian.df <- data.frame(t(apply(vertices, 1, function(row_) barycentric2cartesian(row_, v_cartesian))))

g <- create_plot_triangle(graphsettings) +
  geom_point(data=vertices_cartesian.df, aes(x = X1, y = X2))
g

## and we can plot it as a polygon
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=vertices_cartesian.df, aes(x = X1, y = X2))
g
## but here we see that the polygon is looking weird, because we haven't cared about the ordering of the points.
## If we care we obtain a better looking polygon:
average_ <- colMeans(vertices_cartesian.df)
o_ <- order(apply(sweep(vertices_cartesian.df, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
vertices_cartesian.df <- vertices_cartesian.df[o_,]
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=vertices_cartesian.df, aes(x = X1, y = X2))
g

## next we check whether randomly drawn points are inside the polytope or not
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=vertices_cartesian.df, aes(x = X1, y = X2))
g
## draw points
unif_samples <- gtools::rdirichlet(1e4, c(1,1,1))
unif_samples_cartesian <- t(apply(unif_samples, 1, function(row) barycentric2cartesian(row, v_cartesian)))
within_ <- apply(unif_samples, 1, function(row) all(eta_linearconstraints(etas)$constr %*% row <= eta_linearconstraints(etas)$rhs))
unif_samples_cartesian.df <- data.frame(X1 = unif_samples_cartesian[1:1000,1], X2 = unif_samples_cartesian[1:1000,2], within = within_[1:1000])
g + geom_point(data = unif_samples_cartesian.df, aes(x = X1, y = X2, col = within)) + viridis::scale_color_viridis(discrete=T)

## proportion of points within
mean(within_)

###### 
# library(volesti)
# P = Hpolytope$new(cvxpolytope$constr$constr, cvxpolytope$constr$rhs)
# volume(P, Algo = "CG", error = 0.0001)
## 0.004443476

## no match but maybe because this is not adjusting for the volume of the encompassing simplex
# Psimplex = Hpolytope$new(cvxpolytope$constr$constr[1:3,], cvxpolytope$constr$rhs[1:3])
# volume(Psimplex, Algo = "CG", error = 0.0001)
## very close to 0.5

## so 
# mean(within_)
## should match 
# 2 * volume(P, Algo = "CG", error = 0.0001)



