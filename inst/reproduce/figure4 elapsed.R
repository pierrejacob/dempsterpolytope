## This scripts times the iterations of the Gibbs sampler,
## for various K and N, and plots the timings.

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


Ks <- c(4,8,12,16)
samplesizes <- c(256, 512, 1024, 2048)

## record timings for different K and N 
df_ <- data.frame()
repeats <- 50
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

