## This scrips times some calculations to do with the proposed method.
## It was used to write the rejoinder of the article.

rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)
theme_set(ggthemes::theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), 
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1), 
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1), 
             legend.text = element_text(size = 20), 
             legend.title = element_text(size = 20), title = element_text(size = 30), 
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"), 
             legend.position = "bottom")

K <- 50
N <- 500
counts <- tabulate(sample(x=1:K, size = N, replace = TRUE), nbins = K)
print(sum(counts>0))

lag <- 30
nrep <- 100
meeting_times <- unlist(foreach(irep = 1:nrep) %dorng% sample_meeting_times(counts, lag = lag))
tmax <- floor(max(meeting_times)*1.2)
ubounds <- sapply(1:tmax, function(t) tv_upper_bound(meeting_times, lag, t))
plot(x = 1:tmax, y = ubounds, type = 'l', xlab = "iteration", ylab = "TV upper bounds")

mixing_time_TV1pct <- which(ubounds < 0.01)[1]
cat("time to get to stationarity (1% in TV):", mixing_time_TV1pct, "MCMC iterations\n")

niterations <- 10
repeats <- 10
elapseds <- c()
for (irepeat in 1:repeats){
  pct <- proc.time()
  gibbs_results <- gibbs_sampler(niterations = niterations, counts = counts)
  elapseds <- c(elapseds, (proc.time() - pct)[3])
}
median(elapseds)
cat("time to perform", niterations, "MCMC iterations:", median(elapseds), "seconds (median over", repeats, "runs)\n")

## add empty categories
niterations <- 150
gibbs_results <- gibbs_sampler(niterations, counts)
extended_counts <- c(counts, rep(0, 150))
extended_K <- length(extended_counts)
repeats <- 10
elapseds <- c()
for (irepeat in 1:repeats){
  pct <- proc.time()
  extended_gibbs_results <- extend_us(gibbs_results$Us, whichbefore = 1:K, whichnew = (K+1):extended_K)
  elapseds <- c(elapseds, (proc.time() - pct)[3])
}
print(elapseds)
print(median(elapseds))
cat("time to extend number of categories:", median(elapseds), "seconds (median over", repeats, "runs)\n")

extended_gibbs_results <- extend_us(gibbs_results$Us, whichbefore = 1:K, whichnew = (K+1):extended_K)
dim(extended_gibbs_results$etas)
dim(extended_gibbs_results$etas[1,,])

burnin <- 50
niterations <- dim(extended_gibbs_results$etas)[1]
etas <- extended_gibbs_results$etas[(burnin+1):niterations,,]

## Solve linear programs
extended_K <- dim(etas)[2]
objvec <- rep(0, extended_K)
objvec[1] <- 1
repeats <- dim(etas)[1]
elapseds <- c()
min1st <- c()
for (irepeat in 1:repeats){
  pct <- proc.time()
  min1st <- c(min1st, lpsolve_over_eta(etas[irepeat,,], objvec))
  elapseds <- c(elapseds, (proc.time() - pct)[3])
}
print(summary(elapseds))
print(median(elapseds))
hist(min1st)
cat("time to solve LP:", median(elapseds), "seconds (median over", repeats, "runs)\n")

## Solve quadratic programs

library(quadprog)
##
mindist <- function(eta, phat){
  K <- dim(eta)[1]
  ## solve QP to find closest point to P in set
  ## we want to minimize (y - p)^T (y - p) = y^T y + p^T p - 2 p^T y
  ## which is equivalent to minimizing y^T y - 2 p^T y 
  ## which is equivalent to minimizing 0.5 y^T D y - d^T y with D = diag(2) and d = 2p
  Dmat <- diag(2, K, K)
  dvec <- 2 * phat
  ## constraints: given by simplex
  
  # number of constraints in the LP: K+1 constraints for the simplex
  # and K*(K-1) constraints of the form theta_ell / theta_k < eta_{k,ell}
  nconstraints <- (K + 1) + K*(K-1)
  # matrix encoding the constraints
  mat_cst <- matrix(0, nrow = nconstraints, ncol = K)
  mat_cst[1,] <- 1
  for (i in 1:K) mat_cst[1+i,i] <- 1
  # direction of constraints
  dir_ <- c("=", rep(">=", K), rep(">=", K*(K-1)))
  # right hand side of constraints
  rhs_ <- c(1, rep(0, K), rep(0, K*(K-1)))
  row <- K+1
  for (k in 1:K){
    for (ell in setdiff(1:K, k)){
      row <- row + 1
      ## constraint of the form
      # theta_i - eta_{j,i} theta_j <= 0 
      if (all(is.finite(eta[k,]))){
        mat_cst[row,k] <- eta[k,ell]
        mat_cst[row,ell] <- -1
      }
    }
  }
  ## solve QP
  solution.QP <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(mat_cst), bvec = rhs_, meq = 1)
  ## compute smallest distance
  l <- sum((phat - solution.QP$solution)^2)
  return(l)
}

phat <- extended_counts / sum(extended_counts)
repeats <- dim(etas)[1]
elapseds <- c()
mindists <- c()
for (irepeat in 1:repeats){
  pct <- proc.time()
  mindists <- c(mindists, mindist(etas[irepeat,,], phat))
  elapseds <- c(elapseds, (proc.time() - pct)[3])
}
print(summary(elapseds))
print(median(elapseds))
hist(mindists)
cat("time to solve QP:", median(elapseds), "seconds (median over", repeats, "runs)\n")

### don't run!
## etas_vertices(etas[1,,])

