rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)

counts <- c(8,4)
nsamples <- 30
exact_dirichlet_pts <- gtools::rdirichlet(nsamples, counts+c(2,2))
cvxpolytope_cartesian.df <- data.frame()
for (isample in 1:nsamples){
  cvx_Vrepresentation <- rbind(c(exact_dirichlet_pts[isample,], 0), c(0,0,1))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, isample = isample))
}
## 
g <- create_plot_triangle(graphsettings = graphsettings)
gdirichlet <- g + geom_path(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = isample), alpha = .25, lineend = 'round')
gdirichlet
##
# ggsave(plot = gdirichlet, filename = "partialprior1.pdf", width = 5, height = 5)

K <- 3
new_counts <- c(2,1,3)
new_niterations <- 1e3
new_results <- gibbs_sampler(new_niterations, new_counts)
subiter <- floor(seq(from = 1e2, to = new_niterations, length.out = 30))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(new_results$etas[iter,,])
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 0.25, alpha = .25, fill = 'black', colour = 'black')
gpolytopes
# gdirichletpolytopes <- gdirichlet + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), alpha = .6, fill = 'black', colour = 'black')

# ggsave(plot = gpolytopes, filename = "partialprior2.pdf", width = 5, height = 5)



