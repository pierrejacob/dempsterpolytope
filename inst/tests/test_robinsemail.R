# devtools::install_github("pierrejacob/dempsterpolytope")
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())
## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
set.seed(4)
##
categories <- 1:K
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
gtriangle <- g

## initially we have two categories
counts <- c(3,2)
niterations <- 1e3
results <- gibbs_sampler(niterations, counts)

##
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 20))
offset <- truncnorm::rtruncnorm(n = niterations, a = -.02, b = .02, sd = .01)
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(results$etas[iter,,])
  cvx_Vrepresentation <- cbind(cvx$vertices_barcoord, 0)
  cvx_Vrepresentation <- rbind(cvx_Vrepresentation,
                               cbind(cvx$vertices_barcoord, offset[iter]))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gintervals <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                                       alpha = .4)
gintervals

cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  etas <- matrix(1, nrow = 3, ncol = 3)
  etas[1:2,1:2] <- results$etas[iter,,]
  etas[1,3] <- etas[2,3] <- Inf
  etas[3,1] <- etas[3,2] <- Inf
  cvx <- etas2cvxpolytope(etas)
  cvx_Vrepresentation <- cvx$vertices_barcoord
  cvx_Vrepresentation[1:2,3] <- offset[iter]
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}

gintervalext <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 2,
                                         alpha = .4)
gintervalext

## what if we run the Gibbs sampler with data (3,2,0)
counts <- c(3,2,0)
niterations <- 1e3
results <- gibbs_sampler(niterations, counts)
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 20))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(results$etas[iter,,])
  cvx_Vrepresentation <- cbind(cvx$vertices_barcoord, 0)
  cvx_Vrepresentation <- rbind(cvx_Vrepresentation,
                               cbind(cvx$vertices_barcoord, offset[iter]))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gpolytopes <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                                       alpha = .4)
gpolytopes

## it's different because knowing N_3 = 0 is information on theta_3 (explicitly in the form of upper bounds on it)


## let's go back to 2 categories
## here we have an explicit form for eta_{k to l} as a ratio of Beta variables

counts <- c(3,2,1)
niterations <- 1e3
results <- gibbs_sampler(niterations, counts[1:2])
dim(results$Us[[1]])
dim(results$Us[[2]])
iter <- 500
results$Us[[1]][iter,,]
results$Us[[2]][iter,,]
## Note: 
## the 1st component of the three Us in category 1 is less than that of two Us in category 2
## that's because these three Us are "opposite" vertex 1 

## we can reconstruct etas 
## with one on diagonal
## and eta_{1,2} = min_{u in category 1} u_{2}/u_{1}
results$etas[iter,1,2]
min(results$Us[[1]][iter,,2] / results$Us[[1]][iter,,1])
## and eta_{2,1} = min_{u in category 2} u_{1}/u_{2}
results$etas[iter,2,1]
min(results$Us[[2]][iter,,1] / results$Us[[2]][iter,,2])


## we can view this as points on the interval
## and the following vector is distributed as sorted uniforms
c(sort(results$Us[[1]][iter,,1]), sort(results$Us[[2]][iter,,1]))

## so if we collect e.g. over the iterations
sortedUs_gibbs <- foreach(irep = 1:niterations, .combine = rbind) %dopar% c(sort(results$Us[[1]][irep,,1]), sort(results$Us[[2]][irep,,1]))
## we know that its marginal distribution is that of a sorted uniform
sortedUs_exact <- apply(matrix(runif(1e4*sum(counts[1:2])), ncol = 1e4), 2, sort)
sortedUs_exact[,1:3]
##
imarg <- 3
hist(sortedUs_gibbs[,imarg], prob = TRUE, nclass = 100)
hist(sortedUs_exact[imarg,], prob = TRUE, add = TRUE, col = rgb(1,0,0,0.5), nclass = 100)
curve(dbeta(x, imarg, sum(counts[1:2])-imarg+1), add=T)
##

## what if we create these intervals on each pair (k,ell)
## "minimally extend" and then take the intersection?

## simulate eta_{k,ell} directly
etas_12s <- foreach(irep = 100:niterations, .combine = c) %dopar% results$etas[irep,1,2]
hist(log(etas_12s), prob = TRUE, nclass = 100)
x <- rbeta(1e4, counts[1], counts[2] + 1)
hist(log((1-x)/x), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)

etas_21s <- foreach(irep = 100:niterations, .combine = c) %dopar% results$etas[irep,2,1]
hist(log(etas_21s), prob = TRUE, nclass = 100)
x <- rbeta(1e4, counts[1]+1, counts[2])
hist(log((x/(1-x))), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)

## uniforms for categories k,l
etas_minextintersection <- function(counts){
  reject <- TRUE
  while(reject){
    etas <- diag(1, 3, 3)
    pairs <- list()
    pairs[[1]] <- c(1,2)
    pairs[[2]] <- c(1,3)
    pairs[[3]] <- c(2,3)
    for (pair in pairs){
      k <- pair[1]
      ell <- pair[2]
      gamma_k <- rgamma(1, counts[k], 1)
      gamma_ell <- rgamma(1, counts[ell], 1)
      exp_kell <- rexp(n = 1, rate = 1)
      etas[k,ell] <- (gamma_ell + exp_kell) / gamma_k
      etas[ell,k] <- (gamma_k + exp_kell) / gamma_ell
    }
    isnonempty  <- dempsterpolytope:::check_cst_graph(etas)
    reject <- !isnonempty
  }
  return(etas)
}

etas <- etas_minextintersection(counts)
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

triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
g <- g + scale_color_discrete("category: ")
for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = contcols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = contcols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}
# g

cvx <- etas2cvxpolytope(etas)
cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
average_ <- colMeans(cvx_cartesian)
o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
cvx_cartesian <- cvx_cartesian[o_,]
g <- g + geom_polygon(data = data.frame(cvx_cartesian), aes(x = X1, y = X2), alpha = .4)
g

netas <- 1e4
etas_minextintersection_results <- array(dim = c(netas, 3, 3))
for (irep in 1:netas){
  etas_minextintersection_results[irep,,] <- etas_minextintersection(counts)
}

### first question, does the rejection mechanism affect the marginal distribution of etas_{kl}?
### we know what it was before the rejection mechanism, so we can compare.
hist(log(etas_minextintersection_results[,1,2]), prob = TRUE, nclass = 100, ylim = c(0,.6))
x <- rbeta(1e4, counts[1], counts[2] + 1)
hist(log((1-x)/x), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)
## it does not seem so!

## another of looking at it is to plot the marginal distribution of etas_{kl}
## when generated by the Gibbs sampler, and compare with what we get if we had only counts k and l
niterations_gibbs <- 1e4
burnin <- 1e3
gibbs_results <- gibbs_sampler(niterations_gibbs, counts)
etas_gibbs <- gibbs_results$etas[(burnin+1):niterations_gibbs,,]

etas_12_gibbs <- etas_gibbs[,1,2]
hist(log(etas_12_gibbs), prob = TRUE, nclass = 100, ylim = c(0,.6))
x <- rbeta(1e5, counts[1], counts[2] + 1)
hist(log((1-x)/x), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)
hist(log(etas_minextintersection_results[,1,2]), prob = TRUE, nclass = 100, add = TRUE, col = rgb(.2,1,.2,0.5))
## the Gibbs sampler produces etas such that etas[1,2] is as if there were only 
## two categories and counts N1 and N2
## and the intersection of the three "minimum extension" does the same thing

## next, compare calculation of some lower/upper probabilities

## marginal lower/upper cdf
grid01 <- seq(from = 0, to = 1, length.out = 50)
imarg <- 3
loweruppercdf_gibbs <- etas_to_lower_upper_cdf_dopar(etas_gibbs, imarg, grid01)
lowercdf_gibbs <- colMeans(loweruppercdf_gibbs$iscontained)
uppercdf_gibbs <- colMeans(loweruppercdf_gibbs$intersects)
loweruppercdf_minext <- etas_to_lower_upper_cdf_dopar(etas_minextintersection_results, imarg, grid01)
lowercdf_minext <- colMeans(loweruppercdf_minext$iscontained)
uppercdf_minext <- colMeans(loweruppercdf_minext$intersects)
# plot lower and upper CDFs
plot(x = grid01, y = lowercdf_gibbs, type = "l", xlab = expression(theta), ylab = "CDF")
lines(x = grid01, y = uppercdf_gibbs)
lines(x = grid01, y = lowercdf_minext, col = "red")
lines(x = grid01, y = uppercdf_minext, col = "red")

## so all indicates that there is a difference between what the Gibbs sampler
## is doing and what the intersections of three min ext are doing

## next what about Dempster's proposal? By design it will be respecting the marginal 
## distribution of etas[k,l], if I understand correctly. Let's check.

counts <- c(3,2,1)
dempster_CRN_sampler <- function(counts){
  ## draw Gamma for each vertex, with as shape the associated count
  gamma_vertices <- sapply(counts, function(shape) rgamma(n = 1, shape = shape, rate = 1))
  ## 
  names(gamma_vertices) <- 1:K
  exp_edges <- c()
  etas <- diag(1, 3, 3)
  for (k in 1:(K-1)){
    for (ell in setdiff(k:K, k)){
      gamma_k <- gamma_vertices[paste0(k)]
      gamma_ell <- gamma_vertices[paste0(ell)]
      exp_ <- rexp(n = 1, rate = 1)
      exp_edges <- c(exp_edges, exp_)
      names(exp_edges)[length(exp_edges)] <- paste0(k, "-", ell)
      # interval_k_to_ell <- c(gamma_k/(gamma_k+exp_+gamma_ell), (gamma_k+exp_)/(gamma_k+exp_+gamma_ell))
      # etaktoell <- interval_k_to_ell[2] / (1 - interval_k_to_ell[2])
      # etaelltok <- (1 - interval_k_to_ell[1]) / interval_k_to_ell[1]
      etas[k,ell] <- (gamma_ell + exp_) / gamma_k
      etas[ell,k] <- (gamma_k + exp_) / gamma_ell
    }
  }
  return(etas)
}

netas <- 1e4
etas_CRN <- array(dim = c(netas, 3, 3))
for (irep in 1:netas){
  etas_CRN[irep,,] <- dempster_CRN_sampler(counts)
}

subiter <- floor(seq(from = 2e3, to = netas, length.out = 50))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(etas_CRN[iter,,])
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gpolytopes <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                                       alpha = .4)
gpolytopes

nsides_CRN <- c()
for (iter in 1:netas){
  cvx <- etas2cvxpolytope(etas_CRN[iter,,])
  nsides_CRN <- c(nsides_CRN,nrow(cvx$vertices_barcoord))
}
table(nsides_CRN)/length(nsides_CRN)

##
nsides_gibbs <- c()
for (iter in (burnin+1):niterations_gibbs){
  cvx <- etas2cvxpolytope(gibbs_results$etas[iter,,])
  nsides_gibbs <- c(nsides_gibbs,nrow(cvx$vertices_barcoord))
}
table(nsides_gibbs)/length(nsides_gibbs)

nsides_minext <- c()
for (iter in (burnin+1):niterations_gibbs){
  cvx <- etas2cvxpolytope(etas_minextintersection_results[iter,,])
  nsides_minext <- c(nsides_minext,nrow(cvx$vertices_barcoord))
}
table(nsides_minext) / length(nsides_minext)


# dempsterpolytope:::check_cst_graph(etas_CRN[10,,])
# hist(log(etas_CRN[,1,2]), prob = TRUE, nclass = 100, ylim = c(0,.6))
# x <- rbeta(1e5, counts[1], counts[2] + 1)
# hist(log((1-x)/x), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)
# 
# hist(log(etas_CRN[,2,1]), prob = TRUE, nclass = 100, ylim = c(0,.6))
# x <- rbeta(1e5, counts[2], counts[1] + 1)
# hist(log((1-x)/x), add=T, col = rgb(1,0,0,0.5), prob=T, nclass=100)

loweruppercdf_CRN <- etas_to_lower_upper_cdf_dopar(etas_CRN, imarg, grid01)
lowercdf_CRN <- colMeans(loweruppercdf_CRN$iscontained)
uppercdf_CRN <- colMeans(loweruppercdf_CRN$intersects)

plot(x = grid01, y = lowercdf_gibbs, type = "l", xlab = expression(theta), ylab = "CDF")
lines(x = grid01, y = uppercdf_gibbs)
lines(x = grid01, y = lowercdf_minext, col = "red")
lines(x = grid01, y = uppercdf_minext, col = "red")
lines(x = grid01, y = lowercdf_CRN, col = "blue")
lines(x = grid01, y = uppercdf_CRN, col = "blue")




