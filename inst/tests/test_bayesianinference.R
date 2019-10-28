library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

# number of observations 
# (don't set it to be too large, because rejection sampler would become very slow; values of <= 5,6 are OK)
n <- 6
# number of categories
K <- 3
categories <- 1:K
## for plotting purposes
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]

# data 
# counts <- c(70,150,0)
# theta_dgp <- c(0.3, 0.3, 0.4)
# X <- sample(x = categories, size = n, replace = TRUE, prob = theta_dgp)
X <- c(categories, sample(x = categories, size = n - K, replace = TRUE))
counts <- tabulate(X)
counts

## we put a Dirichlet prior distribution 
## hyper parameter
alpha <- rep(1/K, K)
## then postior
post_parameters <- alpha + counts
## sample from posterior 
rpost <- function(nsamples){
  gtools::rdirichlet(nsamples, post_parameters)
}

postsamples <- rpost(1e4)
#
par(mfrow = c(1,3))
hist(postsamples[,1])
hist(postsamples[,2])
hist(postsamples[,3])
par(mfrow = c(1,1))
hist(postsamples[,1], prob = TRUE, ylim = c(0, 5), col = rgb(0,0,0))
hist(postsamples[,2], prob = TRUE, col = rgb(1,0,0,0.5), add = TRUE)
hist(postsamples[,3], prob = TRUE, col = rgb(1,1,0,0.5), add = TRUE)

## plot the samples in a triangle
# convert points to cartesian coordinates
postsamples_cartesian <- t(apply(postsamples, 1, function(v) barycentric2cartesian(v, v_cartesian)))
library(ggpointdensity)
library(viridis)
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
gpoints <- g + geom_pointdensity(data = data.frame(x = postsamples_cartesian[,1], y = postsamples_cartesian[,2]), mapping = aes(x = x, y = y)) + 
  scale_color_viridis() + theme(legend.position = "none")
gpoints 


## how to plot the pdf 
# grid on the simplex
# expand grid with two sequences from 0 to 1 
ngrid <- 120
df <- expand.grid(x1 = seq(from = 0, to = 1, length.out = ngrid), x2 = seq(from = 0, to = 1, length.out = ngrid))
df <- df %>% mutate(w1 = x1) %>% mutate(w2 = x2) %>% filter(w1 + w2 <= 1) %>% mutate(w3 = 1-w1-w2) 
df_cartesian <- data.frame(t(apply(df %>% select(starts_with("w")), 1, function(v) barycentric2cartesian(v, v_cartesian))))
# plot(df_cartesian[,1], df_cartesian[,2])
df_cartesian$z <- apply(df, 1, function(v) gtools::ddirichlet(v[3:5], post_parameters))
head(df_cartesian)
g + geom_point(data = df_cartesian, aes(x = X1, y = X2, colour = z, fill = z)) + scale_color_viridis() + theme(legend.position = "none")

## animation showing difference as prior changes

df_cart_z <- data.frame()
for (al in seq(from = c(0.1, to = 10, length(10)))){
  alpha <- rep(al, K)
  post_parameters <- alpha + counts
  z <- apply(df, 1, function(v) gtools::ddirichlet(v[3:5], post_parameters))
  df_cart_z <- rbind(df_cart_z, data.frame(X1 = df_cartesian$X1, X2 = df_cartesian$X2, al = al, z))
}

tail(df_cart_z)
library(gganimate)
ggplot(df_cart_z, aes(x = X1, y = X2, colour = z, fill = z)) + geom_point() + scale_color_viridis() + theme(legend.position = "none") + transition_states(al) + 
  labs(title="state: {closest_state}")

## then postior

##






