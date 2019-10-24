library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]

####
set.seed(4)
# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(9,8,3)

## run Gibbs sampler
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, counts = counts)
iteration <- 100
etas <- samples_gibbs$etas[iteration,,]
pts_barcoord <- lapply(samples_gibbs$Us, function(l) l[iteration,,])


pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))
g <- ggplot_triangle(v_cartesian, pts_cart, etas, addpolytope = T, cols = cols)
g <- g + geom_label(data = data.frame(x = .73, y = 1, label = TeX("$\\theta_1/\\theta_3 = \\eta_{3\\rightarrow 1}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = "blue", alpha = 0.75)
g <- g + geom_label(data = data.frame(x = 0.8, y = 0.78, label = TeX("$\\theta_3/\\theta_1 = \\eta_{1\\rightarrow 3}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = "red", alpha = 0.75)
g <- g + geom_label(data = data.frame(x = 0.15, y = 0.7, label = TeX("$\\theta_2/\\theta_1 = \\eta_{1\\rightarrow 2}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = "red", alpha = 0.75)
g

ggsave(filename = "sdk.plottriangle.points.pdf", plot = g, width = 5, height = 5)


## now get all polytopes of feasible parameters at all iterations
## and overlay them in plot
df.polytope <- data.frame()
for (iteration in 100:niterations){
  etas <- samples_gibbs$etas[iteration,,]
  etascvxp <- etas2cvxpolytope(etas)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  df.polytope <- rbind(df.polytope, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], 
                                               iteration = iteration))
}
g <- ggplot_triangle(v_cartesian) +
  geom_polygon(data=df.polytope %>% filter(iteration >= 100), aes(x = x, y = y, group = iteration), alpha = .3)
g

ggsave(filename = "sdk.plottriangle.polytopeS.pdf", plot = g, width = 5, height = 5)


