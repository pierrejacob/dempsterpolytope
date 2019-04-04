library(montecarlodsm)
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
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(pi[1](theta)), parse = FALSE, size = 5) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(pi[2](theta)), parse = FALSE, size = 5) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(pi[3](theta)), parse = FALSE, size = 5) 
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.pdf", plot = g, width = 5, height = 5)

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
g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[1], linetype = 2)
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
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.linearconstraint.pdf", plot = g, width = 5, height = 5)


####
set.seed(4)
# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
freqX <- c(9,8,3)
niterations <- 200
samples_gibbs <- gibbs_sampler(niterations = niterations, freqX = freqX)
##
iteration <- 100
etas <- samples_gibbs$etas_chain[iteration,,]
pts_barcoord <- lapply(samples_gibbs$Achain, function(l) l[iteration,,])
pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))
g <- ggplot_triangle(v_cartesian, pts_cart, etas, addpolytope = T, cols = cols)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.points.pdf", plot = g, width = 5, height = 5)


####

# function to plot triangle with ggplot
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
for (d in categories){
  g <- g + geom_point(data=data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2]), colour = cols[d])
}
g <- g + scale_color_discrete("category: ")

for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}

## constraint of the form eta_12 eta_21 >= 1 
# theta_1/theta_2 < eta_21
# eta_12^{-1} < theta_1/theta_2 
add_L2const <- function(g, i1, i2, etas){
  A <- matrix(rep(1, K-1), ncol = K-1)
  A <- rbind(A, diag(-1, K-1, K-1))
  b <- c(1, rep(0, K-1))
  ccc <- rep(0, K)
  ccc[i1] <- 1 
  ccc[i2] <- -etas[i2,i1]
  cc <- ccc - ccc[K]
  b <- c(b, -ccc[K])
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  ccc <- rep(0, K)
  ccc[i1] <- -1 
  ccc[i2] <- etas[i1,i2]^{-1}
  cc <- ccc - ccc[K]
  b <- c(b, -ccc[K])
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  vertices_barcoord <- hitandrun::findVertices(constr)
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
  return(g)
}
g <- add_L2const(g, 1, 2, etas)
g <- add_L2const(g, 1, 3, etas)
g <- add_L2const(g, 2, 3, etas)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.constraintsL2.pdf", plot = g, width = 5, height = 5)


# function to plot triangle with ggplot
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
for (d in categories){
  g <- g + geom_point(data=data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2]), colour = cols[d])
}
g <- g + scale_color_discrete("category: ")

for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}
add_L3const <- function(g, i1, i2, i3, etas){
  ## constraint of the form eta_12 eta_23 eta_31 >= 1 
  # theta_1/theta_2 < eta_21 
  # theta_2/theta_3 < eta_32
  # theta_3/theta_1 < eta_13
  A <- matrix(rep(1, K-1), ncol = K-1)
  A <- rbind(A, diag(-1, K-1, K-1))
  b <- c(1, rep(0, K-1))
  ccc <- rep(0, K)
  ccc[i1] <- 1 
  ccc[i2] <- -etas[i2,i1]
  cc <- ccc - ccc[K]
  b <- c(b, -ccc[K])
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  ccc <- rep(0, K)
  ccc[i2] <- 1 
  ccc[i3] <- -etas[i3,i2]
  cc <- ccc - ccc[K]
  b <- c(b, -ccc[K])
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  ccc <- rep(0, K)
  ccc[i3] <- 1 
  ccc[i1] <- -etas[i1,i3]
  cc <- ccc - ccc[K]
  b <- c(b, -ccc[K])
  A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  vertices_barcoord <- hitandrun::findVertices(constr)
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
  return(g)
}

g <- add_L3const(add_L3const(g, 1, 2, 3, etas), 3, 2, 1, etas)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.constraintsL3.pdf", plot = g, width = 5, height = 5)

### Constraint violations

triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
g <- g + scale_color_discrete("category: ")
etas1 <- etas
etas1[2,3] <- 5
etas1[3,1] <- 2

for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas1[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas1[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}
g <- add_L3const(add_L3const(g, 1, 2, 3, etas1), 3, 2, 1, etas1)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.violateconstraintsL2.pdf", plot = g, width = 5, height = 5)



triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
g <- g + scale_color_discrete("category: ")
etas1 <- etas
etas1[1,2] <- 10 
etas1[2,1] <- 0.2 

for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas1[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas1[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}
# g <- add_L3const(add_L3const(g, 1, 2, 3, etas1), 3, 2, 1, etas1)
g <- add_L2const(g, 1, 2, etas1)
g <- add_L2const(g, 1, 3, etas1)
g <- add_L2const(g, 2, 3, etas1)
g
# ggsave(filename = "~/Dropbox/Fiducial/Figures/sdk.plottriangle.violateconstraintsL3.pdf", plot = g, width = 5, height = 5)
