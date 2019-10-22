## This scripts computes upper bounds on the total variation distance
## between the marginal distribution of the Gibbs sampler,
## for various K and N, and plots the upper bounds as a function of the iteration.

library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

#### Effect of the number of categories on the mixing time
## draw meeting times

## increase NREP to get smoother curves in the resulting figures
NREP <- 5e2
# arbitrary choice of lag (see Biswas et al.)
lag <- 75
omega <- 0.9
meetings <- list()
# number of categories
Ks <- c(5, 10, 20)

for (iK in seq_along(Ks)){
  K <- Ks[iK]
  cat("K =", K, "\n")
  freqX <- rep(10, K)
  n <- sum(freqX)
  meetings[[iK]] <- unlist(foreach(irep = 1:NREP) %dorng% {
    meeting_times(freqX, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
  })
}
save(Ks, meetings, file = "inst/reproduce/mixing.K.RData")
load(file = "inst/reproduce/mixing.K.RData")

niterations <- 100
ubounds.df <- data.frame()
for (iK in seq_along(Ks)){
  ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetings[[iK]], lag, t))
  ubounds.df <- rbind(ubounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds, K = Ks[iK]))
}


g <- ggplot(ubounds.df, aes(x = iteration, y = ubounds, group = K, linetype = factor(K))) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g <- g + scale_linetype(name = "K")
g

ggsave(plot = g, filename = "inst/reproduce/mixing.K.pdf", width = 7, height = 5)
# g + scale_y_log10()

#### Effect of the number of observations on the mixing time
###
lag_multiplier <- 5
omega <- 0.9
meetings <- list()
# number of categories
K <- 5
Ns <- c(10, 20, 30, 40)
for (iN in seq_along(Ns)){
  N <- Ns[iN]
  lag <- lag_multiplier * N
  cat("N =", N, "\n")
  freqX <- rep(N, K)
  n <- sum(freqX)
  meetings[[iN]] <- unlist(foreach(irep = 1:NREP) %dorng% {
    meeting_times(freqX, lag = lag, rinit = function(){ x = rexp(K); return(x/sum(x))}, omega = omega, max_iterations = 1e5)
  })
}
save(Ns, lag_multiplier, K, meetings, file = "inst/reproduce/mixing.N.RData")
load(file = "inst/reproduce/mixing.N.RData")

niterations <- 200
ubounds.df <- data.frame()
for (iN in seq_along(Ns)){
  ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetings[[iN]], lag_multiplier * Ns[iN], t))
  ubounds.df <- rbind(ubounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds, N = K*Ns[iN]))
}


g <- ggplot(ubounds.df, aes(x = iteration, y = ubounds, group = N, linetype = factor(N))) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g <- g + scale_linetype(name = "N")
g
ggsave(plot = g, filename = "inst/reproduce/mixing.N.pdf", width = 7, height = 5)

# scaling of mixing time with N 
ubounds.df %>% group_by(N) %>% summarise(tmix = iteration[which(ubounds < 0.01)[1]])
g + geom_hline(yintercept = 0.01, linetype = 1)

