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

K <- 3
counts <- c(1,2,1)
N <- sum(counts)

niterations <- 5e3
burnin <- 9e2
samples_gibbs <- dempsterpolytope::gibbs_sampler(niterations, counts)
names(samples_gibbs)


## now get all polytopes of feasible parameters at all iterations
## and overlay them in plot
df.polytope <- data.frame()
for (iteration in (burnin+1):niterations){
  etas <- samples_gibbs$etas[iteration,,]
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
g <- create_plot_triangle(graphsettings) +
  geom_polygon(data=df.polytope %>% filter(iteration >= 1000, iteration <= 1100), 
               aes(x = x, y = y, group = iteration), alpha = .2, colour = 'black')
g

### Monte Carlo estimation of the volume
## let's take a convex polytope
etas <- samples_gibbs$etas[niterations,,]
etascvxp <- etas2cvxpolytope(etas)
names(etascvxp)


## sample uniformly in the simplex
unif_samples <- gtools::rdirichlet(1e3, c(1,1,1))
head(unif_samples)

library(volesti)
# ?volesti::Vpolytope
# volesti::Vpolytope$new()
volumes <- c()
for (iter in (burnin+1):niterations){
  if (iter %% 100 == 1) print(iter)
  cvx_ <- etas2cvxpolytope(samples_gibbs$etas[iter,,])
  # A = cvx_$constr$constr
  # b = cvx_$constr$rhs
  unif_samples <- gtools::rdirichlet(1e3, c(1,1,1))
  within_ <- apply(unif_samples, 1, function(row) all(cvx_$constr$constr %*% row[1:2] <= cvx_$constr$rhs))
  volumes <- c(volumes, mean(within_))
  # A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
  # b = c(0,0,1)
  # P = Hpolytope$new(A, b)
  # volumes <- c(volumes, 2 * volesti::volume(P, Algo = 'CG', error = 0.1))
}

mean(volumes)
1/choose(N+K-1,K-1)

