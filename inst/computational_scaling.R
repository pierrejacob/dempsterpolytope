library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

samplesizes <- c(500, 1000, 1500, 2000)
Ks <- c(3,5,7,9)
df_ <- data.frame()
for (n in samplesizes){
  for (K in Ks){
    categories <- 1:K
    # data-generating parameter
    theta_dgp <- rexp(K)
    theta_dgp <- theta_dgp / sum(theta_dgp)
    # data 
    X <- sample(x = categories, size = n - K, replace = TRUE, prob = theta_dgp)
    X <- c(X, categories)
    # frequencies
    freqX <- as.numeric(tabulate(X, nbins = K))
    ###
    niterations <- 100
    pct <- proc.time()
    samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
    elapsed <- (proc.time() - pct)[3]
    df_ <- rbind(df_, data.frame(n = n, K = K, elapsed = elapsed))
    ##
  }
}

head(df_)
g <- ggplot(df_, aes(x = n, y = elapsed, group = K, colour = factor(K))) + geom_line() + ggthemes::scale_color_colorblind(name = "K:")
g <- g + theme(legend.position = "bottom") + xlab("N")
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.scalingn.pdf", plot = g, width = 5, height = 5)

g <- ggplot(df_, aes(x = K, y = elapsed, group = n, colour = factor(n))) + geom_line() + ggthemes::scale_color_colorblind(name = "N:")
g <- g + theme(legend.position = "bottom")
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.scalingK.pdf", plot = g, width = 5, height = 5)


## test what happens if we had empty categories
K <- 3
n <- 100
categories <- 1:K
# data-generating parameter
theta_dgp <- rexp(K)
theta_dgp <- theta_dgp / sum(theta_dgp)
# data 
X <- sample(x = categories, size = n - K, replace = TRUE, prob = theta_dgp)
X <- c(X, categories)
# frequencies
freqX <- as.numeric(tabulate(X, nbins = K))
###
niterations <- 100
pct <- proc.time()
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
elapsed <- (proc.time() - pct)[3]
elapsed
freqX_ <- freqX
df_ <- data.frame()
for (addK in 1:50){
  freqX_ <- c(freqX_, 0)
  pct <- proc.time()
  samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX_)
  elapsed <- (proc.time() - pct)[3]
  df_ <- rbind(df_, data.frame(n = n, K = K+addK, elapsed = elapsed))
}

g <- ggplot(df_, aes(x = K, y = elapsed)) + geom_line() 
g <- g + theme(legend.position = "bottom")
g
