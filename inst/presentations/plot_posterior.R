library(dempsterpolytope)
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

counts <- c(2,3,1)
## we put a Dirichlet prior distribution 
## hyper parameter
alpha <- rep(1, K)
## then postior
post_parameters <- alpha + counts
## plot the samples in a triangle
# convert points to cartesian coordinates
library(ggpointdensity)
library(viridis)
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
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
gdensity <- g + geom_point(data = df_cartesian, aes(x = X1, y = X2, colour = z, fill = z)) + scale_color_viridis() + theme(legend.position = "none")
ggsave(filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/presentation/posteriordensity.png", plot = gdensity, width = 5, height = 5)
