library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
library(gridExtra)
set.seed(1)
rm(list = ls())

# number of observations
n <- 50
# number of categories
K <- 3
categories <- 1:K
# data 
# freqX <- c(70,150,0)
# theta_dgp <- c(0.3, 0.3, 0.4)
# X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)

X <- c(categories, sample(x = categories, size = n - K, replace = TRUE))
freqX <- tabulate(X)
freqX

## encompassing triangle has three vertices
if (K == 3){
  v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
  cols <- c("red", "green", "blue")
  v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
  matrixT <- matrix(0, nrow = K-1, ncol = K-1)
  for (k in 1:(K-1)){
    kth_components <- sapply(v_cartesian, function(x) x[k])
    matrixT[k,] <- kth_components[-K] - kth_components[K]
  }
  #
  ## show constraints
  barconstraint2cartconstraint <- function(d, j, eta, matrixT, v_cartesian){
    # ccc * (wA wB wC) = 0 with:
    ccc <- rep(0, 3)
    ccc[d] <- 1
    ccc[j] <- - eta
    # which is equivalent to ftilde * (wA wB) = gtilde with
    ftilde <- c(ccc[1] - ccc[3], ccc[2] - ccc[3])
    gtilde <- -ccc[3]
    # we can generically express that as a constraint on x,y through
    f <- solve(t(matrixT), ftilde)
    g <- gtilde + sum(f * v_cartesian[[3]])
    # f1 x + f2 y = g is equivalent to a = g/f2, b = - f1/f2
    return(c(g/f[2], -f[1]/f[2]))
  }
}



###
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)

##
length(samples_gibbs$Achain)
dim(samples_gibbs$Achain[[1]])
dim(samples_gibbs$Achain[[2]])
dim(samples_gibbs$Achain[[3]])

if (K == 3){
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
}

### now SMC sampler
nparticles <- 2^5

set.seed(1)
smc_res_graph <- SMC_sampler_graph(nparticles, X, K, verbose = TRUE)
set.seed(1)
smc_res_lp <- SMC_sampler_lp(nparticles, X, K, verbose = TRUE)

smc_res_graph$etas_particles[1,,]
smc_res_lp$etas_particles[1,,]


# etas <- smc_res$etas_particles[1,,]
# g <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
# g2 <- g
# E(g2, c(1,2))$weight <- 0.5
# E(g2, c(1,3))$weight <- 0.5
# 
# E(g2, c(1,2))$weight
# E(g2, c(1,3))$weight
# 
# g3 <- g
# E(g3, c(1,2,1,3))$weight <- c(0.5, 0.5)
# E(g3, c(1,2))$weight
# E(g3, c(1,3))$weight
# k_ <- 2
# minimum_values <- rep(1, K)
# for (ell in setdiff(1:K, k_)){
#   minimum_values[ell] <- distances(g, v = ell, to = k_, mode = "out")
# }
# minimum_values
# minimum_values2 <- rep(1, K)
# minimum_values2[setdiff(1:K, k_)] <-distances(g, v = setdiff(1:K, k_), to = k_, mode = "out")[,1]


library(microbenchmark)
microbenchmark(
  smcgraph = SMC_sampler_graph(nparticles, X, K),
  smclp = SMC_sampler_lp(nparticles, X, K),
  times = 10)


nrep <- 2*(detectCores()-2)
smc_res <- foreach(irep = 1:nrep) %dorng% {
  SMC_sampler_graph(nparticles, X, K)
}

var(sapply(smc_res, function(v) sum(v$normcst)))
hist(sapply(smc_res, function(v) sum(v$normcst)))

smc_res_lp <- foreach(irep = 1:nrep) %dorng% {
  SMC_sampler_lp(nparticles, X, K)
}

var(sapply(smc_res_lp, function(v) sum(v$normcst)))
hist(sapply(smc_res_lp, function(v) sum(v$normcst)))



## now try to assimilate X[2]
# assimilate <- function(etas, x){
#   k_ <- x

#  visualize etas
plot_etas <- function(etas){
  triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
  g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
  for (d in categories){
    # set indices for two other components
    j1 <- setdiff(categories, d)[1]
    j2 <- setdiff(categories, d)[2]
    interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
    interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
    g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
    g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
    intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
  }
  g
}
if (K==3){
  grid.arrange(plot_etas(smc_res_graph$etas_particles[1,,]), plot_etas(smc_res_graph$etas_particles[2,,]),
             plot_etas(smc_res_graph$etas_particles[3,,]), plot_etas(smc_res_graph$etas_particles[4,,]), nrow = 2)
}


smc_res <- SMC_sampler_lp(2^10, X, K, verbose = FALSE)
etas_particles <- smc_res$etas_particles
weights <- smc_res$weights

### test
burnin <- 1000
niterations <- 5000 + burnin
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
param <- 1
interval <- c(freqX[1]/sum(freqX)-0.05, freqX[1]/sum(freqX)+0.05)
intervalcvxp <- interval2polytope(K, param, interval)
# 
postburn <- niterations - burnin
contained_mcmc <- rep(0, postburn)
intersects_mcmc <- rep(0, postburn)
for (index in ((burnin+1):niterations)){
  cvxp <- etas2cvxpolytope(samples_gibbs$etas_chain[index,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_mcmc[index-burnin] <- res_[1]
  intersects_mcmc[index-burnin] <- res_[2]
}

contained_smc <- rep(0, length(weights))
intersect_smc <- rep(0, length(weights))
for (iparticle in 1:length(weights)){
  cvxp <- etas2cvxpolytope(etas_particles[iparticle,,])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_smc[iparticle] <- weights[iparticle] * res_[1]
  intersect_smc[iparticle] <- weights[iparticle] * res_[2]
}


# for the interval
interval
# on parameter
param
# equal to 
(freqX/sum(freqX))[param]
# lower/upper probabilities, post burn-in
cat(mean(contained_mcmc), mean(intersects_mcmc), "\n")
cat(sum(contained_smc), sum(intersect_smc), "\n")

