## This scripts times the iterations of the Gibbs sampler,
## for various K and N, and plots the timings.

library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# samplesizes <- c(256, 512)
# Ks <- c(4,8)
samplesizes <- c(256, 512, 1024, 2048)
Ks <- c(4,8,12,16)

df_ <- data.frame()
repeats <- 10
for (n in samplesizes){
  for (K in Ks){
    counts <- rep(floor(n/K), K)
    niterations <- 100
    elapseds <- c()
    for (irepeat in 1:repeats){
      pct <- proc.time()
      samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
      elapseds <- c(elapseds, (proc.time() - pct)[3])
    }
    elapsed <- median(elapseds)
    df_ <- rbind(df_, data.frame(n = n, K = K, elapsed = elapsed))
  }
}

save(df_, file = "computational.scaling.RData")
load(file = "computational.scaling.RData")

g <- ggplot(df_, aes(x = n, y = elapsed, group = K, linetype = factor(K))) + geom_line() 
g <- g + scale_linetype(name = "K") + xlab("N") + scale_x_continuous(breaks = samplesizes)
g
ggsave(filename = "elapsed.scalingn.pdf", plot = g, width = 7, height = 5)

g <- ggplot(df_, aes(x = K, y = elapsed, group = n, linetype = factor(n))) + geom_line() 
g <- g + scale_linetype(name = "N") + scale_x_continuous(breaks = Ks)
g
ggsave(filename = "elapsed.scalingK.pdf", plot = g, width = 7, height = 5)

