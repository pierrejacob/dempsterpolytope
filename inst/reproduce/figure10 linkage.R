## This scripts reproduce the results for the linkage model

rm(list = ls())
set.seed(2)
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
theme_set(ggthemes::theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), 
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1), 
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1), 
             legend.text = element_text(size = 20), 
             legend.title = element_text(size = 20), title = element_text(size = 30), 
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"), 
             legend.position = "bottom")


# define A and b such that theta = A phi + b where phi is in the interval (0,1)
A <- c(1/4, -1/4, -1/4, 1/4)
b <- c(1/2, 1/4, 1/4, 0)

# data 
K <- 4
counts = c(25, 3, 4, 7)

## 1. Dirichlet DSM model.

# function to sample n times from Dirichlet(alpha) distribution
# (taken verbatim from the 'gtools' package)
rdirichlet <- function (n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

## function to generate samples from the Dirichlet DSM (DDSM) model 
sample_ddsm <- function(x, nsim = 1e4, batch = nsim*100){
  # Rejection algorithm for Dirichlet DSM
  # x: 4 dimensional data vector
  # nsim: total samples wanted
  # batch: sample size per batch
  
  sample_phi <- function(n, data = x){
    z = rdirichlet(n, alpha = c(data, 1))
    s_ = cbind(4*pmax(z[, 1]-1/2, z[, 4]), 1 - 4*pmax(z[, 2], z[, 3]))
    return(s_[s_[,2]>s_[,1],])
  }
  
  s_ <- sample_phi(batch, data = x); ns = nrow(s_)
  if (ns == 0){
    stop('Acceptance rate too low, change batch size')
  }
  
  while(ns < nsim){
    s_ <- rbind(s_, sample_phi(batch, data = x))
    ns <- nrow(s_)
  }
  return(s_[1:nsim,])
}

## draw 1e4 feasible sets
phi_intervals_ddsm = sample_ddsm(x = counts, nsim = 1e4)
## show empirical lower/upper CDF
plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red',
     main = 'Dirichlet DSM upper and lower CDF, nsim = 1e4',
     xlab = 'phi0', ylab = 'P(phi < phi0)', xlim = c(0, 1))
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)

## 2. Imprecise Dirchlet Model.
post_idm <- function(x, a, by = 1e-3){
  # trapezoid integration for idm posterior
  # x: data; a: alpha
  # function returns a vector of cdf evaluated per ``by''
  phi <- seq(0, 1, by = by)
  p_ <- (1/2+phi/4)^(x[1]+a[1]-1)*
    (1/4-phi/4)^(x[2]+x[3]+a[2]+a[3]-1)*
    (phi/4)^(x[4]+a[4]-1)
  vol_ <- (p_[-1] + head(p_, -1))*by/2
  cdf <- cumsum(vol_)/sum(vol_)
  return(c(0, cdf))
}

# a grid of alpha values to try
a_grid <- expand.grid(a1 = seq(0, 1, by = 0.02),
                      a2 = seq(0, 1, by = 0.02),
                      a4 = seq(0, 1, by = 0.02))
a_grid <- a_grid[which(rowSums(a_grid) == 1), ]
a_grid$a3 <- 0
a_grid <- as.matrix(a_grid[, c(1,2,4,3)])

# evaluate idm posterior at a grid of a values
by = 1e-4
idm_grid <- array(NA, dim = c(nrow(a_grid), 1/by+1))
for (i in 1:nrow(a_grid)){
  alpha <- a_grid[i, ]
  idm_grid[i, ] <- post_idm(x=counts, a=alpha, by=by)
}

# pointwise lower and upper probabilities of IDM
lb_idm <- apply(idm_grid, 2, min)
ub_idm <- apply(idm_grid, 2, max)

plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'DDSM vs IDM upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
lines(seq(0, 1, by=by), lb_idm, col = 'blue', lty = 2)
lines(seq(0, 1, by=by), ub_idm, col = 'red', lty = 2)
legend('topleft', lty = c(1, 2), legend = c('DDSM', 'IDM'))

## 3. Simplex DSM model.
omega <- 0.9
lag <- 100
NREP <- 1e3
meetingtimes <- foreach(irep = 1:NREP, .combine = c) %dorng% {
  meeting_times(counts, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
}
## TV upper bounds
niterations <- 100
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes, lag, t))
g <- ggplot(data = data.frame(iteration = 1:niterations, ubounds = ubounds), aes(x = iteration, y = ubounds)) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g 
## let's choose a burn-in of 50, based on the above plot
burnin <- 50
nchains <- 250
niterations <- 500
library(abind)
acomb <- function(...) abind(..., along=1)
etas <- foreach(irep = 1:nchains, .combine = 'acomb') %dorng% {
  init <- rexp(K)
  init <- init/sum(init)
  samples_gibbs <- gibbs_sampler(niterations, counts, theta_0 = init)
  samples_gibbs$etas[(burnin+1):niterations,,]
}
dim(etas)
##
## define segment constraints in dimension K-1
constr <- matrix(0, nrow = 2, ncol = K-1)
rhs <- c(-b[1], A[1] + b[1])
constr[1,1] <- -1
constr[2,1] <- 1
dir <- c("<=", "<=")
for (j in 2:(K-1)){
  row_constr_ <- rep(0, K-1)
  row_constr_[1] <- -A[j]/A[1]
  row_constr_[j] <- 1
  rhs <- c(rhs, b[j] - b[1] * A[j] / A[1], -(b[j] - b[1] * A[j] / A[1]))
  constr <- rbind(constr, row_constr_, -row_constr_, deparse.level = 0)
  dir <- c(dir, "<=", "<=")
}
segment_constr <- list(constr = constr, dir = dir, rhs = rhs)

intersects_with_segment <- function(eta, segment_constr){
  cvxp <- etas2cvxpolytope(eta)
  test_constr <- segment_constr
  test_constr$constr <- rbind(test_constr$constr, cvxp$constr$constr)
  test_constr$rhs <- c(test_constr$rhs, cvxp$constr$rhs)
  test_constr$dir <- c(test_constr$dir, cvxp$constr$dir)
  ## make H representation (H for ?)
  h <- rcdd::makeH(test_constr$constr, test_constr$rhs)
  ## try to find V representation (V for Vendetta or Vertex?)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  intersects <- (dim(v)[1] != 0)
  phis_ <- c(NA, NA)
  if (intersects){
    seg <- v[,-c(1,2)]
    seg <- cbind(seg, 1-rowSums(seg))
    phis_ <- apply(seg, 1, function(v) (v[1] - b[1]) / A[1])
    phis_ <- (sort(phis_))
  }
  phis_
}

intersections_ <- foreach(iter = 1:dim(etas)[1], .combine = rbind) %dopar% {
  intersects_with_segment(etas[iter,,], segment_constr)
}

dim(intersections_)
## acceptance rate
mean(!is.na(intersections_[,1]))
## retain feasible sets that intersect with segment constraint 
intersections_ <- intersections_[!is.na(intersections_[,1]),]

linkage_cdf_lowerupper <- function(phi_0, lu_chain){
  # intersect if phi_lower < phi_0
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < phi_0))
  # contains if phi_upper < phi_0
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < phi_0))
  return(c(cdf_lower, cdf_upper))
}

phi_grid_length <- 50
phi_grid <- seq(from = 0, to = 1, length.out = phi_grid_length)
cdf_lower_values <- rep(0, phi_grid_length)
cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], intersections_)
  cdf_lower_values[iphi] <- cdf_res[1]
  cdf_upper_values[iphi] <- cdf_res[2]
}
##

plot(ecdf(phi_intervals_ddsm[, 1]), col = 'red', xlim = c(0, 1),
     main = 'DDSM vs IDM vs SDSM upper and lower CDF',
     xlab = 'phi0', ylab = 'P(phi < phi0)')
plot(ecdf(phi_intervals_ddsm[, 2]), col = 'blue', add = T)
lines(seq(0, 1, by=by), lb_idm, col = 'blue', lty = 2)
lines(seq(0, 1, by=by), ub_idm, col = 'red', lty = 2)
lines(phi_grid, cdf_lower_values, col = 'blue', lty = 3)
lines(phi_grid, cdf_upper_values, col = 'red', lty = 3)
legend('topleft', lty = c(1, 2, 3), legend = c('DDSM', 'IDM', 'SDSM'))

nddsm <- nrow(phi_intervals_ddsm)
nsdsm <- nrow(intersections_)

df_ <- data.frame(x = c(phi_intervals_ddsm[,1], phi_intervals_ddsm[,2], intersections_[,1], intersections_[,2]),
                  side = c(rep("lower", nddsm), rep("upper", nddsm), rep("lower", nsdsm), rep("upper", nsdsm)),
                  dsm = c(rep("Dirichlet", 2*nddsm), rep("Simplex", 2*nsdsm)))


g <- ggplot(df_, aes(x = x, linetype = dsm, group = interaction(side, dsm))) + stat_ecdf()
g <- g + xlab(expression(phi)) +  ylab("cdf") + scale_linetype(name = "DSM: ")
g

ggsave(filename = "linkage.cdf.pdf", plot = g, width = 6, height = 4)

##
ddsm_cdf_lower_values <- rep(0, phi_grid_length)
ddsm_cdf_upper_values <- rep(0, phi_grid_length)
for (iphi in 1:phi_grid_length){
  ddsm_cdf_res <- linkage_cdf_lowerupper(phi_grid[iphi], phi_intervals_ddsm)
  ddsm_cdf_lower_values[iphi] <- ddsm_cdf_res[1]
  ddsm_cdf_upper_values[iphi] <- ddsm_cdf_res[2]
}
plot(phi_grid, cdf_upper_values - cdf_lower_values, type = "l", ylim = c(0,.17), lty = 2)
lines(phi_grid, ddsm_cdf_upper_values - ddsm_cdf_lower_values)
gr <- qplot(x = phi_grid, y =  ddsm_cdf_upper_values - ddsm_cdf_lower_values, geom = "line", linetype = "Dirichlet") + 
  geom_line(aes(y = cdf_upper_values - cdf_lower_values, linetype = "Simplex")) + ylab("r-cdf") + xlab(expression(phi))
gr <- gr + scale_linetype(name = "DSM: ")
gr            

ggsave(filename = "linkage.r.pdf", plot = gr, width = 6, height = 4)



