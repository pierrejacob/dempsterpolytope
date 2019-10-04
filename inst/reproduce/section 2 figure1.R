## This script generates the two plots of Figure 1 of the article.
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
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
##
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
#

pt_bar <- c(0.3, 0.25, 0.45)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)
# Now with ggplot2
g <- ggplot_triangle(v_cartesian)
g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v2[1], v3[1]), y = c(pt_xy[2], v2[2], v3[2])),
                      aes(fill = "1"), colour = "black", alpha = 0.5)
g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v3[1]), y = c(pt_xy[2], v1[2], v3[2])),
                      aes(fill = "2"), colour = "black", alpha = 0.5)
g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v2[1]), y = c(pt_xy[2], v1[2], v2[2])),
                      aes(fill = "3"), colour = "black", alpha = 0.5)
g <- g + scale_fill_discrete(name = "partition: ") + theme(legend.position = "none")
meanpi1 <- c(pt_xy[1]/3 + v2[1]/3 + v3[1]/3, pt_xy[2]/3 + v2[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(Delta[1](theta)), parse = FALSE, size = 5) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(Delta[2](theta)), parse = FALSE, size = 5) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(Delta[3](theta)), parse = FALSE, size = 5) 
g <- g + annotate(geom = 'text', size = 5, x = .56, y = .28, label = TeX("$\\theta$", output = "character"), parse = TRUE)
g
# ggsave(filename = "sdk.plottriangle.pdf", plot = g, width = 5, height = 5)

## show constraints
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

g <- ggplot_triangle(v_cartesian)
# g <- g + geom_point(aes(x = pt_xy[1], y = pt_xy[2]))
interslope_j1 <- barconstraint2cartconstraint(1, 3, pt_bar[1]/pt_bar[3], matrixT, v_cartesian)
g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[3], linetype = 2)
# g <- g + annotate(geom = "text", x = pt_xy[1]-0.025, y = pt_xy[2]+0.025, label = expression(theta), parse = FALSE, size = 5) 

## add polytope theta_3/theta_1 < pt_bar[3]/pt_bar[1]
## ie - pt_bar[3]/pt_bar[1] theta_1 + theta_3 < 0
A <- matrix(rep(1, K-1), ncol = K-1)
A <- rbind(A, diag(-1, K-1, K-1))
b <- c(1, rep(0, K-1))
# then we add the constraint theta_k < b
ccc <- rep(0, K)
ccc[1] <- pt_bar[3]/pt_bar[1]
ccc[3] <- -1
cc <- ccc - ccc[K]
b <- c(b, -ccc[K])
A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))

constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
vertices_barcoord <- hitandrun::findVertices(constr)
vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
# g <- g + annotate(geom = 'text', size = 5, x = .9, y = .48, label = TeX("$\\theta_1/\\theta_3 = \\eta$", output = "character"), parse = TRUE)
g <- g + annotate(geom = 'text', size = 5, x = .7, y = .18, label = TeX("$\\theta_1/\\theta_3 \\leq \\eta$", output = "character"), parse = TRUE)
g

# ggsave(filename = "sdk.plottriangle.linearconstraint.pdf", plot = g, width = 5, height = 5)




