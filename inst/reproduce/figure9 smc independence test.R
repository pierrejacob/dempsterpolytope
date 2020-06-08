### This scripts reproduces the examples of Section 7 of the paper on independence testing
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


counts <- c(16, 5, 14, 18)
K <- 4
counts_table = matrix(counts, byrow = T, nrow = 2)

posterior_association <- function(data, prior = 1, nsim = 1e5){
  r_ = mapply(function(l){rgamma(nsim, l)}, data+prior)
  positive = mean(r_[,1]*r_[,4] > r_[,2]*r_[,3])
  return(data.frame(positive = positive, negative = 1-positive, 
                    se = sqrt(positive*(1-positive)/nsim)))
}
posterior_association(data = counts_table)


## test function
# h <- function(etas) unlist(compare_with_independence(etas))
nchains <- 100
niterations <- 200
burnin <- 100

ds_estimators <- foreach(irep = 1:nchains, .combine = rbind) %dorng% {
  init <- rexp(K)
  init <- init/sum(init)
  samples_gibbs <- gibbs_sampler(niterations, counts, theta_0 = init)
  etas <- samples_gibbs$etas[(burnin+1):niterations,,]
  positivassoc <- t(apply(etas, 1, positiveassociation))
  c(mean(positivassoc[,1]), mean(positivassoc[,2]))
}

positivassoc <- colMeans(ds_estimators)
positivassoc.p <- positivassoc[2]
positivassoc.r <-  positivassoc[1] - positivassoc[2]
positivassoc.q <- 1 - positivassoc.p - positivassoc.r
cat(positivassoc.p, positivassoc.q, positivassoc.r)


## SMC approach 
X <- unlist(sapply(seq_along(counts), function(n) rep(n, counts[n])))
X <- sample(x = X, size = length(X), replace = F)
X

nparticles <- 2^9
samples_smc <- SMC_sampler(nparticles, X, K, essthreshold = 0.75, verbose = TRUE)
names(samples_smc)


smc_pqr <- foreach (iobs = 1:length(X), .combine = rbind) %dopar% {
  positivassoc <- t(apply(samples_smc$etas_history[[iobs]], 1, positiveassociation))
  positivassoc <- colSums(positivassoc * samples_smc$weights_history[[iobs]])
  positivassoc.p <- positivassoc[2]
  positivassoc.r <-  positivassoc[1] - positivassoc[2]
  positivassoc.q <- 1 - positivassoc.p - positivassoc.r
  c(positivassoc.p, positivassoc.q, positivassoc.r)
}

head(smc_pqr)
# matplot(cbind(smc_pqr[,1], smc_pqr[,1]+smc_pqr[,3]), type = 'l')

## create data frame to add to background
pqr.df <- data.frame(x = seq_along(X), ymin = smc_pqr[,1], ymax = smc_pqr[,1]+smc_pqr[,3], row.names = NULL)
background.df <- data.frame(row.names = NULL)
for (i in 2:length(X)){
  if (X[i] %in% c(1,4)){
    background.df <- rbind(background.df, data.frame(i = i, x = c(i-1, i, i, i-1), y = c(0, 0, smc_pqr[i,1], smc_pqr[i-1,1]), col = "14"))
  } else {
    # background.df <- rbind(background.df, data.frame(xmin = i-1, xmax = i, ymin = 0, ymax = smc_pqr[i,1], col = "23"))
  }
}
row.names(background.df) <- NULL
g <- ggplot(data = pqr.df,
            aes(x = x, ymin = ymin, ymax = ymax))
g <- g + geom_polygon(data = background.df, aes(x = x, y = y, group = factor(i), ymin = NULL, ymax = NULL), fill = 'grey')
g <- g + geom_ribbon() + scale_x_continuous(breaks = seq_along(X), labels = X)
g <- g + xlab('observations') + ylab('positive association')
g <- g + scale_fill_manual(values = c('grey', 'white')) + theme(legend.position = 'bottom')
g
ggsave(filename = "sequentialpositiveassociation.pdf", plot = g, width = 14, height = 4)

