library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations
n <- 3
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
##
## encompassing triangle has three vertices
if (K == 3){
  v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
  cols <- c("red", "green", "blue")
}

###
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
##
samples_gibbs$etas_chain[100,,]

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

### Rejection sampling
# sample a_n uniformly in simplex
sample_etas <- function(){
  a <- matrix(rexp(K*n), ncol = K)
  a <- t(apply(a, 1, function(v) v / sum(v)))
  etas <- diag(1, K, K)
  for (k in 1:K){
    notk <- setdiff(1:K, k)
    a_k <- a[X == k,,drop=F]
    for (j in notk){
      etas[k,j] <- min(a_k[,j]/a_k[,k])
    }
  }
  etas
}
# check constraints
# with graph tools
check_cst_graph <- function(etas){
  g <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
  negativecycle <- inherits(try(distances(g, mode = "out"), silent = TRUE), "try-error")
  return(!negativecycle)
}

# or manually, for K = 3
# check_cst_manual <- function(etas){
#   satis <- TRUE
#   for (k in 1:K){
#     notk <- setdiff(1:K, k)
#     for (j in notk){
#       l <- setdiff(1:K, c(k,j))
#       if (etas[k,j] < 1/(etas[j,k])){
#         satis <- FALSE
#       }
#       if (etas[k,j] < 1/(etas[j,l] * etas[l,k])){
#         satis <- FALSE
#       }
#     }
#   }
#   return(satis)
# }

NREP <- 1e4
etas_ <- foreach(irep = 1:NREP) %dorng% {
  etas <- sample_etas()
  sat <- check_cst_graph(etas)
  list(etas = etas, sat = sat)
}

mean(sapply(etas_, function(x) x$sat))

## check with SMC
nparticles <- 2^10
smc_res <- SMC_sampler(nparticles, X, K)
# log normalizing constant
exp(sum(smc_res$normcst))

etas_accepted <- etas_[which(sapply(etas_, function(x) x$sat))]
if (K == 3){
  rejection.df.polytope <- data.frame()
  for (iteration in 1:length(etas_accepted)){
    etas <- etas_accepted[[iteration]]$etas
    etascvxp <- etas2cvxpolytope(etas)
    ## convert coordinates to cartesian
    vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
    # order vertices according to angles
    average_ <- colMeans(vertices_cart)
    o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
    vertices_cart <- vertices_cart[o_,]
    rejection.df.polytope <- rbind(rejection.df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                                                     iteration = iteration))
  }
  g <- ggplot_triangle(v_cartesian) +
    geom_polygon(data=rejection.df.polytope, aes(x = x, y = y, group = iteration), alpha = .3)
  g
}  
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.polytopeS.pdf", plot = g, width = 7, height = 7)


##
niterations <- 5100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
elapsed <- (proc.time() - pct)[3]
elapsed
burnin <- 100

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

## now with rejection sampling
nrs <- 5e3
etas_rs <- foreach(irep = 1:nrs) %dorng% {
  accept <- FALSE
  etas <- NULL
  while (!accept){
    etas <- sample_etas()
    accept <- check_cst_graph(etas)
  }
  etas
}
contained_rs <- rep(0, nrs)
intersects_rs <- rep(0, nrs)
for (index in 1:nrs){
  cvxp <- etas2cvxpolytope(etas_rs[[index]])
  res_ <- compare_polytopes(cvxp, intervalcvxp)
  contained_rs[index] <- res_[1]
  intersects_rs[index] <- res_[2]
}

# lower/upper probabilities, obtained with rejection sampling
cat(mean(contained_rs), mean(intersects_rs), "\n")
# versus the Gibbs output 
cat(mean(contained_), mean(intersects_), "\n")

# seems to agree
