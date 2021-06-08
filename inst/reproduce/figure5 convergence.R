## This scripts computes upper bounds on the total variation distance
## between the marginal distribution of the Gibbs sampler,
## for various K and N, and plots the upper bounds as a function of the iteration.

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


#### Effect of the number of categories on the mixing time
## draw meeting times

## increase NREP to get smoother curves in the resulting figures
NREP <- 5e2
# choice of lag (see Biswas et al.)
lag <- 50
omega <- 0.9
meetings <- list()
# number of categories
# Ks <- c(5, 10, 20)
Ks <- c(4,8,12,16)
for (iK in seq_along(Ks)){
  K <- Ks[iK]
  cat("K =", K, "\n")
  counts <- rep(10, K)
  n <- sum(counts)
  meetings[[iK]] <- unlist(foreach(irep = 1:NREP) %dorng% {
    sample_meeting_times(counts, lag = lag, omega = omega, max_iterations = 1e5)
  })
}
save(Ks, meetings, file = "mixing.K.RData")
load(file = "mixing.K.RData")

niterations <- 100
ubounds.df <- data.frame()
for (iK in seq_along(Ks)){
  ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetings[[iK]], lag, t))
  ubounds.df <- rbind(ubounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds, K = Ks[iK]))
}

g <- ggplot(ubounds.df, aes(x = iteration, y = ubounds, group = K, linetype = factor(K))) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g <- g + scale_linetype(name = "K")
g

ggsave(plot = g, filename = "mixing.K.pdf", width = 6, height = 4)

#### Effect of the number of observations on the mixing time
###
lag_multiplier <- 10
omega <- 0.9
meetings <- list()
# number of categories
K <- 5
Ns <- c(10, 20, 30, 40)
for (iN in seq_along(Ns)){
  N <- Ns[iN]
  lag <- lag_multiplier * N
  cat("N =", N, "\n")
  counts <- rep(N, K)
  n <- sum(counts)
  meetings[[iN]] <- unlist(foreach(irep = 1:NREP) %dorng% {
    sample_meeting_times(counts, lag = lag, omega = omega, max_iterations = 1e5)
  })
}
save(Ns, lag_multiplier, K, meetings, file = "mixing.N.RData")
load(file = "mixing.N.RData")

niterations <- 250
ubounds.df <- data.frame()
for (iN in seq_along(Ns)){
  ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetings[[iN]], lag_multiplier * Ns[iN], t))
  ubounds.df <- rbind(ubounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds, N = K*Ns[iN]))
}

g <- ggplot(ubounds.df, aes(x = iteration, y = ubounds, group = N, linetype = factor(N))) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
g <- g + scale_linetype(name = "N")
g
ggsave(plot = g, filename = "mixing.N.pdf", width = 6, height = 4)

# # scaling of mixing time with N 
# ubounds.df %>% group_by(N) %>% summarise(tmix = iteration[which(ubounds < 0.01)[1]])
# 
