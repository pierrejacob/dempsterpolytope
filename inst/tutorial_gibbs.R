### TO BE DONE: This could become a little tutorial (in rmarkdown?) on how to use the package in a simple example

library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
# freqX <- c(70,150,0)
theta_dgp <- c(0.2, 0.4, 0.4)
X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
freqX <- tabulate(X)
freqX
##
## encompassing triangle has three vertices
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")

###
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX, theta_0 = c(0.8,0.1,0.1))
##

## now get all polytopes of feasible parameters at all iterations
## and overlay them in plot
df.polytope <- data.frame()
for (iteration in 1:niterations){
  etas <- samples_gibbs$etas_chain[iteration,,]
  etascvxp <- etas2cvxpolytope(etas)
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
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.polytopeS.pdf", plot = g, width = 7, height = 7)


### Show an instance of a feasible plot with all six linear constraints
iteration <- 100
etas <- samples_gibbs$etas_chain[iteration,,]
pts_barcoord <- lapply(samples_gibbs$Achain, function(l) l[iteration,,])
pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))
g <- ggplot_triangle(v_cartesian, pts_cart, etas, addpolytope = T, cols = cols)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.polytope.pdf", plot = g, width = 7, height = 7)

##
niterations <- 1100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
elapsed <- (proc.time() - pct)[3]
elapsed
burnin <- 1

param <- 1
interval <- c(0.1, 0.2)
intervalcvxp <- interval2polytope(K, param, interval)
# 
postburn <- niterations - burnin
contained_ <- rep(0, postburn)
intersects_ <- rep(0, postburn)
for (index in ((burnin+1):niterations)){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_[index-burnin] <- res_[1]
  intersects_[index-burnin] <- res_[2]
}

# lower/upper probabilities, post burn-in
cat(mean(contained_), mean(intersects_), "\n")
# for the interval
interval
# on parameter
param
# equal to 
(freqX/sum(freqX))[param]

matplot(cbind(cumsum(contained_)/(1:postburn),cumsum(intersects_)/(1:postburn)), type = "l")
matplot(intersects_)
