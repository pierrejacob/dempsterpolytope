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
set.seed(13)
## number of categories
K <- 3
categories <- 1:K
##
## data 
counts <- c(3,1,2)
## run Gibbs sampler
niterations <- 70
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
(proc.time() - pct)[3]
##

## get a KxK matrix of etas
etas <- samples_gibbs$etas[niterations,,]
etas
## transform that into convex polytope
## where theta is in the simplex and 
## theta_ell / theta_k <= eta[k,l] for all k,l in [K]
cvxpolytope <- etas2vertices(etas)

## vertices of the polytope, in barycentric coordinates
cvxpolytope$vertices_barcoord
## that is, each row represents a vertex, and sums to one

## to plot vertices, we convert the barycentric coordinates to cartesian coordinates
# ?barycentric2cartesian
cvxpolytope_cartesian.df <- data.frame(t(apply(cvxpolytope$vertices_barcoord, 1, function(row_) barycentric2cartesian(row_, v_cartesian))))

g <- create_plot_triangle(graphsettings) +
  geom_point(data=cvxpolytope_cartesian.df, aes(x = X1, y = X2))
g

## and we can plot it as a polygon
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=cvxpolytope_cartesian.df, aes(x = X1, y = X2))
g
## but here we see that the polygon is looking weird, because we haven't cared about the ordering of the points.
## If we care we obtain a better looking polygon:
average_ <- colMeans(cvxpolytope_cartesian.df)
o_ <- order(apply(sweep(cvxpolytope_cartesian.df, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
cvxpolytope_cartesian.df <- cvxpolytope_cartesian.df[o_,]
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=cvxpolytope_cartesian.df, aes(x = X1, y = X2))
g

## next we check whether randomly drawn points are inside the polytope or not
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=cvxpolytope_cartesian.df, aes(x = X1, y = X2))
g
## draw points
unif_samples <- gtools::rdirichlet(1e6, c(1,1,1))
unif_samples_cartesian <- t(apply(unif_samples, 1, function(row) barycentric2cartesian(row, v_cartesian)))
within_ <- apply(unif_samples, 1, function(row) all(cvxpolytope$constr$constr %*% row[1:2] <= cvxpolytope$constr$rhs))
unif_samples_cartesian.df <- data.frame(X1 = unif_samples_cartesian[1:1000,1], X2 = unif_samples_cartesian[1:1000,2], within = within_[1:1000])
g + geom_point(data = unif_samples_cartesian.df, aes(x = X1, y = X2, col = within)) + viridis::scale_color_viridis(discrete=T)

mean(within_)
library(volesti)
P = Hpolytope$new(cvxpolytope$constr$constr, cvxpolytope$constr$rhs)
volume(P, Algo = "CG", error = 0.0001)
## 0.004443476

## no match but maybe because this is not adjusting for the volume of the encompassing simplex
Psimplex = Hpolytope$new(cvxpolytope$constr$constr[1:3,], cvxpolytope$constr$rhs[1:3])
volume(Psimplex, Algo = "CG", error = 0.0001)
## very close to 0.5

## so 
mean(within_)
## should match 
2 * volume(P, Algo = "CG", error = 0.0001)



