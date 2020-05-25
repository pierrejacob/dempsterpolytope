library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_custom_theme()
set.seed(1)
rm(list = ls())

## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
##
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
#

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


####
set.seed(4)
# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(9,8,3)

etas <- structure(c(1, 3.97840569581908, 5.78277269927406, 0.50153243920328, 
                    1, 2.40371666744944, 0.254985206727502, 1.10872229394467, 1), .Dim = c(3L, 
                                                                                           3L))
pts_barcoord <- list(structure(c(0.220719993744219, 0.0814138939835954, 0.182431333296683, 
                                 0.513229103657028, 0.283639301028938, 0.032507041200313, 0.116298407197962, 
                                 0.193913089196759, 0.25781903643718, 0.215218628154736, 0.644720797757343, 
                                 0.357468671333566, 0.355905067248414, 0.179896421174025, 0.448712704448105, 
                                 0.0689457633834079, 0.0972537046182938, 0.215093947968671, 0.564061378101045, 
                                 0.273865308259061, 0.460099995369751, 0.130865829094558, 0.536464277797037, 
                                 0.518780254351582, 0.81475582941863, 0.708833206184947, 0.52708701559415
), .Dim = c(9L, 3L)), structure(c(0.461151703148454, 0.75244214565423, 
                                  0.284243747161101, 0.590621712352687, 0.100832392674065, 0.284007898342033, 
                                  0.590022846892913, 0.349010733144474, 0.0734298522521787, 0.117397086878935, 
                                  0.0714466469469953, 0.0644618442103327, 0.0124837780987313, 0.0182229442420094, 
                                  0.000167534128178599, 0.0537214224488938, 0.465418444599368, 
                                  0.130160767466834, 0.644309605891904, 0.344916443436981, 0.886683829227204, 
                                  0.697769157415957, 0.409809618978908, 0.597267844406632), .Dim = c(8L, 
                                                                                                     3L)), structure(c(0.586172773532202, 0.995516750346673, 0.218482239363358, 
                                                                                                                       0.312461877349755, 0.00316608665436151, 0.753003674024612, 0.101365349118043, 
                                                                                                                       0.00131716299896569, 0.0285140866120302), .Dim = c(3L, 3L)))


pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))


####

# function to plot triangle with ggplot
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
# for (d in categories){
#   g <- g + geom_point(data=data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2]), colour = cols[d])
# }
# g <- g + scale_color_discrete("category: ")

for (d in categories){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = contcols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = contcols[d], linetype = 2)
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
df <- data.frame()
for (d in categories){
  df <- rbind(df, data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2], category = d))
}
df

gpoints <- g + geom_point(data=df, aes(x = x, y = y, fill = factor(category), col = factor(category)), shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpoints

# ggsave(filename = "constraintsL2.pdf", plot = gpoints, width = 5, height = 5)


# # function to plot triangle with ggplot
# triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
# g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
# g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
# g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
# g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
# for (d in categories){
#   g <- g + geom_point(data=data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2]), colour = cols[d])
# }
# g <- g + scale_color_discrete("category: ")
# 
# for (d in categories){
#   # set indices for two other components
#   j1 <- setdiff(categories, d)[1]
#   j2 <- setdiff(categories, d)[2]
#   interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
#   interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
#   g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
#   g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
#   intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
# }
# add_L3const <- function(g, i1, i2, i3, etas){
#   ## constraint of the form eta_12 eta_23 eta_31 >= 1 
#   # theta_1/theta_2 < eta_21 
#   # theta_2/theta_3 < eta_32
#   # theta_3/theta_1 < eta_13
#   A <- matrix(rep(1, K-1), ncol = K-1)
#   A <- rbind(A, diag(-1, K-1, K-1))
#   b <- c(1, rep(0, K-1))
#   ccc <- rep(0, K)
#   ccc[i1] <- 1 
#   ccc[i2] <- -etas[i2,i1]
#   cc <- ccc - ccc[K]
#   b <- c(b, -ccc[K])
#   A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
#   ccc <- rep(0, K)
#   ccc[i2] <- 1 
#   ccc[i3] <- -etas[i3,i2]
#   cc <- ccc - ccc[K]
#   b <- c(b, -ccc[K])
#   A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
#   ccc <- rep(0, K)
#   ccc[i3] <- 1 
#   ccc[i1] <- -etas[i1,i3]
#   cc <- ccc - ccc[K]
#   b <- c(b, -ccc[K])
#   A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
#   constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
#   vertices_barcoord <- hitandrun::findVertices(constr)
#   vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
#   vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
#   g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
#   return(g)
# }
# 
# g <- add_L3const(add_L3const(g, 1, 2, 3, etas), 3, 2, 1, etas)
# g
# 
# ggsave(filename = "constraintsL3.pdf", plot = g, width = 5, height = 5)

