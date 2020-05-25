library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
rm(list = ls())
graphsettings <- set_custom_theme()
set.seed(1)

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
pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, graphsettings$v_cartesian))))

g <- plot_triangle(graphsettings)

barconstraint2cartconstraint <- function(d, j, eta, matrixT, v_cartesian){
  # ccc * (wA wB wC) = 0 with:
  ccc <- rep(0, 3)
  ccc[d] <- 1
  ccc[j] <- - eta
  # which is equivalent to ftilde * (wA wB) = gtilde with
  ftilde <- c(ccc[1] - ccc[3], ccc[2] - ccc[3])
  gtilde <- -ccc[3]
  # we can generically express that as a constraint on x,y through
  f <- solve(t(matrixT), ftilde)
  g <- gtilde + sum(f * v_cartesian[[3]])
  # f1 x + f2 y = g is equivalent to a = g/f2, b = - f1/f2
  return(c(g/f[2], -f[1]/f[2]))
}
for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
  # if (!(removelines)){
    g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = contcols[d], linetype = 2)
    g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = contcols[d], linetype = 2)
  # }
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
  g <- g + geom_polygon(data = data.frame(x = c(v_cartesian[[j1]][1], v_cartesian[[j2]][1], intersection_12[1]),
                                          y = c(v_cartesian[[j1]][2], v_cartesian[[j2]][2], intersection_12[2])), alpha = 0.5, fill = cols[d])
}
for (d in categories){
  g <- g + geom_point(data=data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2]), colour = contcols[d], shape = 21, fill = cols[d])
}
etascvxp <- etas2cvxpolytope(etas)
## convert coordinates to cartesian
vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
# order vertices according to angles
average_ <- colMeans(vertices_cart)
o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
vertices_cart <- vertices_cart[o_,]
g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]))

# g <- g + geom_label(data = data.frame(x = .73, y = 1, label = TeX("$\\theta_1/\\theta_3 = \\eta_{3\\rightarrow 1}$", output = "character")), 
#                     aes(x = x, y = y, label = label), parse = TRUE, col = "blue", alpha = 0.75)
g <- g + geom_label(data = data.frame(x = 0.8, y = 0.78, label = TeX("$\\theta_3/\\theta_1 = \\eta_{1\\rightarrow 3}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = "red", alpha = 0.75)
g <- g + geom_label(data = data.frame(x = 0.15, y = 0.7, label = TeX("$\\theta_2/\\theta_1 = \\eta_{1\\rightarrow 2}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = "red", alpha = 0.75)
g

# ggsave(filename = "gibbsiteration.pdf", plot = g, width = 5, height = 5)


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

# ggsave(filename = "overlaidpolytopes.pdf", plot = g, width = 5, height = 5)


