### This scripts illustrates that inference obtained with or without empty categories are not the same.
### It creates the figures of the rejoinder of the article.

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


# number of categories
K <- 2
categories <- 1:K
# data 
counts <- c(4,3)
cat("Data:", counts, "\n")
##

NREP <- 1000
lag <- 30
meetings <- unlist(foreach(irep = 1:NREP) %dorng% {
  sample_meeting_times(counts, lag = lag)
})

ubounds <- sapply(1:(1.2*lag), function(t) tv_upper_bound(meetings, lag, t))
gtvbounds <- ggplot(data=data.frame(t = seq_along(ubounds), ubounds = ubounds), 
                    aes(x = t, y = ubounds)) + geom_line() + ylab("TV upper bounds") + xlab("iteration")
gtvbounds


niterations <- 1000
burnin <- 100
gibbs_K2 <- gibbs_sampler(niterations, counts)
etas_K2  <- gibbs_K2$etas[(burnin+1):niterations,,]
## find smallest and largest 1st component in polytopes
etas_K2_min1st <- foreach(ieta = 1:(dim(etas_K2)[1]), .combine = c) %dopar% {
  lpsolve_over_eta(etas_K2[ieta,,], c(1,0))
}
etas_K2_max1st <- foreach(ieta = 1:(dim(etas_K2)[1]), .combine = c) %dopar% {
  -lpsolve_over_eta(etas_K2[ieta,,], c(-1,0))
}
##
ecdf_lower_K2 <- ecdf(etas_K2_min1st)
ecdf_upper_K2 <- ecdf(etas_K2_max1st)
##
gfun <- ggplot() + xlim(0,1) + ylim(0,1) + xlab(expression(theta[1])) + ylab("CDF") +
  stat_function(fun = ecdf_lower_K2, colour = 'red', linetype = 1) + stat_function(fun = ecdf_upper_K2, colour = 'red', linetype = 1) 
gfun

### 'r' don't know proba increases
### clean graph
xgrid <- seq(from = 0, to = 1, length.out = 500)
ylower_K2 <- ecdf_lower_K2(xgrid)
yupper_K2 <- ecdf_upper_K2(xgrid)
gribbon <- ggplot() + xlim(0,1) + ylim(0,1) + xlab(expression(theta[1])) + ylab("cdf") + 
  geom_ribbon(data = data.frame(x = xgrid, ymin = ylower_K2, ymax = yupper_K2), aes(x = x, ymin = ymin, ymax = ymax, y = NULL, fill = '2'), alpha = 0.8)
gribbon <- gribbon + scale_colour_manual(name = "K: ", values = c("black", "grey")) + scale_fill_manual(name = "K: ", values = c("black", "grey"))
gribbon

## add more empty categories
df_ <- data.frame(x = xgrid, ymin = ylower_K2, ymax = yupper_K2, K = 2)
library(tictoc)
for (newcat in c(8-2, 32-2, 128-2, 512-2)){
  tic(msg = paste0("K = ", 2 + newcat, " - extension"))
  extended_ <- extend_us(gibbs_K2$Us, whichbefore = c(1,2), whichnew = 2+c(1:newcat))
  extended_etas_ <- extended_$etas
  toc()
  tic(msg = paste0("K = ", 2 + newcat, " - LP"))
  etas_min1st <- foreach(ieta = 1:(dim(extended_etas_)[1]), .combine = c) %dopar% {
    lpsolve_over_eta(extended_etas_[ieta,,], c(1, rep(0, 1 + newcat)))
  }
  etas_max1st <- foreach(ieta = 1:(dim(extended_etas_)[1]), .combine = c) %dopar% {
    -lpsolve_over_eta(extended_etas_[ieta,,], c(-1, rep(0, 1 + newcat)))
  }
  toc()
  ecdf_lower_ <- ecdf(etas_min1st)
  ecdf_upper_ <- ecdf(etas_max1st)
  df_ <- rbind(df_, data.frame(x = xgrid, ymin = ecdf_lower_(xgrid), ymax = ecdf_upper_(xgrid), K = 2+newcat))
}

df_2 <- arrange(df_, desc(K))
df_2$K <- factor(paste0(df_2$K), levels = paste0(unique(df_2$K)),  ordered = TRUE)
df_2$K
gmanyK <- ggplot(data =  df_2, aes(x = x, ymin = ymin, ymax = ymax, 
                                   fill = K, group = K)) +
  geom_ribbon(alpha = 1)
gmanyK <- gmanyK + scale_fill_manual("K:", values = rev(grey.colors(length(unique(df_2$K)))), labels = paste0(unique(df_2$K)),
                                     guide = guide_legend(reverse = TRUE))
gmanyK <- gmanyK + xlim(0,1) + ylim(0,1) + xlab(expression(theta[1])) + ylab("cdf")
gmanyK
ggsave(filename = "rejoinder.emptycategories.theta1cdf.pdf", plot = gmanyK, width = 6, height = 4)

### reproduce the experiments on lower/upper CDF for log (theta_1 / theta_2)

minmax_logtheta1overtheta2 <- function(eta){
  K <- dim(eta)[1]
  results <- rep(0, 2)
  ## find minimum log theta_1 - log theta_2
  vec_ <- rep(0, K)
  vec_[1] <- 1
  vec_[2] <- -1
  if (K > 2){
    vec_[3:K] <- 0
  }
  results[1] <- lpsolve_over_eta_log(eta, vec_)
  results[2] <- -lpsolve_over_eta_log(eta, -vec_)
  return(results)
}


minmax <- foreach(ieta = 1:dim(etas_K2)[1], .combine = rbind) %dopar% {
  minmax_logtheta1overtheta2(etas_K2[ieta,,])
}

ecdf_lower_K2_logtheta1overtheta2 <- ecdf(minmax[,1])
ecdf_upper_K2_logtheta1overtheta2 <- ecdf(minmax[,2])

# lower / upper 
g <- ggplot() + xlim(-3,3) + ylim(0,1) + xlab(expression(log(theta[1]/theta[2]))) + ylab("CDF") +
  stat_function(fun = ecdf_lower_K2_logtheta1overtheta2, colour = 'red', linetype = 1) +
  stat_function(fun = ecdf_upper_K2_logtheta1overtheta2, colour = 'red', linetype = 1)
g

newcat <- 50
extended_ <- extend_us(gibbs_K2$Us, whichbefore = c(1,2), whichnew = 2+c(1:newcat))
extended_etas_ <- extended_$etas

minmax_extended <- foreach(ieta = 1:dim(extended_etas_)[1], .combine = rbind) %dopar% {
  minmax_logtheta1overtheta2(extended_etas_[ieta,,])
}

ecdf_lower_K2_logtheta1overtheta2_extended <- ecdf(minmax_extended[,1])
ecdf_upper_K2_logtheta1overtheta2_extended <- ecdf(minmax_extended[,2])

# lower / upper 
g + stat_function(fun = ecdf_lower_K2_logtheta1overtheta2_extended, colour = 'blue', linetype = 2) +
  stat_function(fun = ecdf_upper_K2_logtheta1overtheta2_extended, colour = 'blue', linetype = 2)

### clean graph
xgrid <- seq(from = -4, to = 4, length.out = 500)
ylower = ecdf_lower_K2_logtheta1overtheta2(xgrid)
yupper = ecdf_upper_K2_logtheta1overtheta2(xgrid)

gratio <- ggplot() + xlim(-4,4) + ylim(0,1) + xlab(expression(log(theta[1]/theta[2]))) + ylab("cdf") + 
  geom_ribbon(data = data.frame(x = xgrid, ymin = ylower, ymax = yupper), aes(x = x, ymin = ymin, ymax = ymax, y = NULL, fill = '2 or 10^69',
                                                                              colour = '2 or 10^69'), alpha = 0.5)
gratio <- gratio  + theme(legend.position = 'bottom') + scale_fill_manual(name = "K: ", values = c("black")) +
  scale_colour_manual(name = "K: ", values = c("black", "grey"))
gratio

ggsave(plot = gratio, filename = "rejoinder.emptycategories.ratio12cdf.pdf", width = 6, height = 4)


