rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)

# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(2,3,1)
# number of observations
n <- sum(counts)

## run Gibbs sampler
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
iteration <- 100
etas <- samples_gibbs$etas[iteration,,]

subiter <- 101:150
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2vertices(samples_gibbs$etas[iter,,])
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
g <- create_plot_triangle(graphsettings = graphsettings)
gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 0.25, alpha = .2, fill = 'black', colour = 'black')
gpolytopes

ggsave(filename = "overlaidpolytopes.pdf", plot = gpolytopes, width = 5, height = 5)


