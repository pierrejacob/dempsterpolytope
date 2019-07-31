library(montecarlodsm)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

samplesizes <- c(256, 512, 1024, 2048)
Ks <- c(4,8,12,16)
Ns <- 
df_ <- data.frame()
repeats <- 10
# for (n in samplesizes){
#   for (K in Ks){
#     freqX <- rep(floor(n/K), K)
#     niterations <- 100
#     elapseds <- c()
#     for (irepeat in 1:repeats){
#       pct <- proc.time()
#       samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
#       elapseds <- c(elapseds, (proc.time() - pct)[3])
#     }
#     elapsed <- median(elapseds)
#     df_ <- rbind(df_, data.frame(n = n, K = K, elapsed = elapsed))
#   }
# }

# save(df_, file = "inst/computational.scaling.RData")
load(file = "inst/computational.scaling.RData")
head(df_)

g <- ggplot(df_, aes(x = n, y = elapsed, group = K, linetype = factor(K))) + geom_line() 
g <- g + scale_linetype(name = "K") + xlab("N")
g
ggsave(filename = "sdk.scalingn.pdf", plot = g, width = 7, height = 5)

g <- ggplot(df_, aes(x = K, y = elapsed, group = n, linetype = factor(n))) + geom_line() 
g <- g + scale_linetype(name = "N")
g
ggsave(filename = "sdk.scalingK.pdf", plot = g, width = 7, height = 5)



# 
# ## test what happens if we had empty categories
# K <- 3
# n <- 100
# categories <- 1:K
# # data-generating parameter
# theta_dgp <- rexp(K)
# theta_dgp <- theta_dgp / sum(theta_dgp)
# # data 
# X <- sample(x = categories, size = n - K, replace = TRUE, prob = theta_dgp)
# X <- c(X, categories)
# # frequencies
# freqX <- as.numeric(tabulate(X, nbins = K))
# ###
# niterations <- 100
# pct <- proc.time()
# samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
# elapsed <- (proc.time() - pct)[3]
# elapsed
# freqX_ <- freqX
# df_ <- data.frame()
# for (addK in 1:50){
#   freqX_ <- c(freqX_, 0)
#   pct <- proc.time()
#   samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX_)
#   elapsed <- (proc.time() - pct)[3]
#   df_ <- rbind(df_, data.frame(n = n, K = K+addK, elapsed = elapsed))
# }
# 
# g <- ggplot(df_, aes(x = K, y = elapsed)) + geom_line() 
# g <- g + theme(legend.position = "bottom")
# g
