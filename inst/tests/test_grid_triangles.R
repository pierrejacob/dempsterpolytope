library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())

## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))

####
set.seed(4)
##
K <- 3
categories <- 1:K
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}

g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
gtriangle <- g
##

counts <- c(1,3,4)
niterations <- 5e3
results <- gibbs_sampler(niterations, counts)
pt_bar <- counts/sum(counts)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)

## function to obtain constraints of the form theta_k / theta_l <= etas_{kl}
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
## function to represent "etas"
plot_etas <- function(etas){
  g <- gtriangle
  for (d in categories){
    # set indices for two other components
    j1 <- setdiff(categories, d)[1]
    j2 <- setdiff(categories, d)[2]
    interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
    interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
    g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = contcols[d], linetype = 2)
    g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = contcols[d], linetype = 2)
    intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
    g <- g + geom_polygon(data = data.frame(x = c(v_cartesian[[j1]][1], v_cartesian[[j2]][1], intersection_12[1]),
                                            y = c(v_cartesian[[j1]][2], v_cartesian[[j2]][2], intersection_12[2])), alpha = 0.2, fill = cols[d])
  }
  etascvxp <- etas2cvxpolytope(etas)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]))
  g <- g + geom_point(data=data.frame(x = pt_xy[1], y = pt_xy[2]), aes(x = x, y = y), colour = 'black', shape = 0, size = 4)
  g <- g + geom_point(data=data.frame(x = pt_xy[1], y = pt_xy[2]), aes(x = x, y = y), colour = 'white', shape = 4, size = 4)
  g
}  

# pl <- lapply(floor(seq(from = 1e3, to = 3e3, length.out = 16)), function(x) plot_etas(results$etas[x,,]))
# ggsave(plot = gridExtra::grid.arrange(grobs = pl, nrow = 4, ncol = 4),
#        filename = "~/Documents/grid1.pdf", width = 30, height = 30)

# pl <- lapply(floor(seq(from = 3e3+1, to = 5e3, length.out = 16)), function(x) plot_etas(results$etas[x,,]))
# ggsave(plot = gridExtra::grid.arrange(grobs = pl, nrow = 4, ncol = 4),
#        filename = "~/Documents/grid2.pdf", width = 30, height = 30)

# pl

plot_etas <- function(etas){
  g <- gtriangle

  etascvxp <- etas2cvxpolytope(etas)
  ## convert coordinates to cartesian
  vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  # order vertices according to angles
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.2)
  g <- g + geom_point(data=data.frame(x = pt_xy[1], y = pt_xy[2]), aes(x = x, y = y), colour = 'black', shape = 0, size = 4)
  g <- g + geom_point(data=data.frame(x = pt_xy[1], y = pt_xy[2]), aes(x = x, y = y), colour = 'white', shape = 4, size = 4)
  
  for (d in categories){
    # set indices for two other components
    j1 <- setdiff(categories, d)[1]
    j2 <- setdiff(categories, d)[2]
    interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
    interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
    intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
    g <- g + geom_polygon(data = data.frame(x = c(v_cartesian[[j1]][1], v_cartesian[[j2]][1], intersection_12[1]),
                                            y = c(v_cartesian[[j1]][2], v_cartesian[[j2]][2], intersection_12[2])), alpha = 0.2, fill = cols[d])
  }
  for (d in categories){
    min_ <- etascvxp$vertices_barcoord[which.min(etascvxp$vertices_barcoord[,d]),]
    max_ <- etascvxp$vertices_barcoord[which.max(etascvxp$vertices_barcoord[,d]),]
    min_cart <- barycentric2cartesian(min_, v_cartesian)
    max_cart <- barycentric2cartesian(max_, v_cartesian)
    g <- g + geom_segment(data = data.frame(x = min_cart[1], xend = max_cart[1], y = min_cart[2], yend = max_cart[2]), 
                          aes(x=x,xend=xend,y=y,yend=yend), colour = cols[d], linetype =d)
  }
  g
  
}  

pl <- lapply(floor(seq(from = 1e3, to = 3e3, length.out = 16)), function(x) plot_etas(results$etas[x,,]))
ggsave(plot = gridExtra::grid.arrange(grobs = pl, nrow = 4, ncol = 4),
       filename = "~/Documents/grid3.pdf", width = 30, height = 30)

pl <- lapply(floor(seq(from = 3e3+1, to = 5e3, length.out = 16)), function(x) plot_etas(results$etas[x,,]))
ggsave(plot = gridExtra::grid.arrange(grobs = pl, nrow = 4, ncol = 4),
       filename = "~/Documents/grid4.pdf", width = 30, height = 30)

# pl

## get the six min/max



