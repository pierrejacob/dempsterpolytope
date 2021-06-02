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
##
K <- 3
categories <- 1:K
gtriangle <- create_plot_triangle(graphsettings)
##
counts <- c(3,2,1)
niterations <- 1e4
results <- gibbs_sampler(niterations, counts)

## now combine this with a Dirichlet prior 
alpha_prior <- c(2,2,2)
prior_pts <- gtools::rdirichlet(niterations, alpha = alpha_prior)
## retain points which fall inside feasible regions
within <- rep(FALSE, niterations)
for (iteration in 1:niterations){
  cvx <- etas2vertices(results$etas[iteration,,])
  within[iteration] <- all(cvx$constr$constr %*% prior_pts[iteration,1:2] <= cvx$constr$rhs)
  # pt_xy <- barycentric2cartesian(prior_pts[i,], v_cartesian)
  # ggplot_triangle(v_cartesian, etas = results$etas[i,,], addpolytope = T) +
  #   geom_point(aes(x = pt_xy[1], y = pt_xy[2]), colour = 'red')
}

mean(within)

post_pts <- prior_pts[within,]
post_pts_cartesian <- t(apply(post_pts, 1, function(row) barycentric2cartesian(row, v_cartesian)))
post_pts_cartesian <- data.frame(post_pts_cartesian)  
gtriangle + geom_point(data = post_pts_cartesian, aes(x = X1, y = X2))

exact_dirichlet_pts <- gtools::rdirichlet(nrow(post_pts), counts+alpha_prior)

par(mfrow = c(1,3))
hist(exact_dirichlet_pts[,1], prob = TRUE, nclass = 40, xlim = c(0,1))
hist(post_pts[,1], prob = TRUE, add = TRUE, nclass = 40, col = rgb(0.6, 0.6, 0.3, 0.5))

hist(exact_dirichlet_pts[,2], prob = TRUE, nclass = 40, xlim = c(0,1))
hist(post_pts[,2], prob = TRUE, add = TRUE, nclass = 40, col = rgb(0.6, 0.6, 0.3, 0.5))

hist(exact_dirichlet_pts[,3], prob = TRUE, nclass = 40, xlim = c(0,1))
hist(post_pts[,3], prob = TRUE, add = TRUE, nclass = 40, col = rgb(0.6, 0.6, 0.3, 0.5))

# qplot(x = exact_dirichlet_pts[,1], geom = "blank") + geom_den


### for illustration purposes
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 100))
cvxpolytope_cartesian.df <- data.frame()
for (iteration in subiter){
  cvxpolytope <- etas2vertices(results$etas[iteration,,])
  cvxpolytope_cartesian <- data.frame(t(apply(cvxpolytope$vertices_barcoord, 1, function(row_) barycentric2cartesian(row_, v_cartesian))))
  average_ <- colMeans(cvxpolytope_cartesian)
  o_ <- order(apply(sweep(cvxpolytope_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvxpolytope_cartesian <- cvxpolytope_cartesian[o_,]
  cvxpolytope_cartesian$iteration <- iteration
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, cvxpolytope_cartesian)
}

gpolytopes <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iteration), alpha = .2, 
                                       colour = 'black')
gpolytopes

prior_pts_cartesian <- t(apply(prior_pts, 1, function(row) barycentric2cartesian(row, v_cartesian)))

gprior <- gtriangle + geom_point(data=data.frame(prior_pts_cartesian[1:1e3,]), aes(x = X1, y = X2))
gprior

gpost <- gtriangle + geom_point(data=data.frame(post_pts_cartesian), aes(x = X1, y = X2))
gpost



