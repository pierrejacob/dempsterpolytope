## This script ...
rm(list = ls())
library(dempsterpolytope)
library(latex2exp)
graphsettings <- set_custom_theme()
set.seed(1)

set.seed(4)
# number of observations
n <- 20
# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(2,3,1)

X <- rep(1:length(counts), times = counts)
K <- length(counts)

sample_uniform <- function(X, K){
  # number of observations
  n <- length(X)
  # matrix of uniform a's in the simplex
  us <- matrix(rexp(K*n), ncol = K)
  us <- t(apply(us, 1, function(v) v / sum(v)))
  pts <- list()
  # create eta
  eta <- diag(1, K, K)
  for (k in 1:K){
    notk <- setdiff(1:K, k)
    u_k <- us[X == k,,drop=F]
    pts[[k]] <- u_k
    for (ell in notk){
      eta[k,ell] <- min(u_k[,ell]/u_k[,k])
    }
  }
  return(list(pts = pts, eta = eta))
}

ntrials <- 1e3
etas <- list()
pts <- list()
accept <- c()
for (itrial in 1:ntrials){
  result <- sample_uniform(X, K)
  etas[[itrial]] <- result$eta
  pts[[itrial]] <- result$pts
  accept <- c(accept, dempsterpolytope:::check_cst_graph(etas[[itrial]]))
}
mean(accept)
index <- which(accept)[1]
eta <- etas[[index]]
pts_barcoord <- pts[[index]]

## draw some points uniformly on the simplex
g <- create_plot_triangle(graphsettings)
gpoints <- add_plot_points(graphsettings, g, pts_barcoord[[1]])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[2]])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[3]])
gpoints 
# ggsave(filename = "stepbystep.feasibleset.1.pdf", plot = gpoints, width = 5, height = 5)

# ## points are each related to one of the x
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[1]], colour = graphsettings$contcols[1], fill = graphsettings$cols[1])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[2]], colour = graphsettings$contcols[2], fill = graphsettings$cols[2])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[3]], colour = graphsettings$contcols[3], fill = graphsettings$cols[3])
gpoints 

eta


# ggsave(filename = "stepbystep.feasibleset.2.pdf", plot = gpoints, width = 5, height = 5)

## show feasible region for theta 
attach(graphsettings)
Asimplex <- matrix(rep(1, K-1), ncol = K-1)
Asimplex <- rbind(Asimplex, diag(-1, K-1, K-1))
bsimplex <- c(1, rep(0, K-1))

polygon_df <- data.frame()

for (k in 1:K){
  A <- Asimplex
  b <- bsimplex
  
  for (othercategory in setdiff(1:K, k)){
    # then we add the constraint theta_k < b
    ccc <- rep(0, K)
    ccc[k] <- -eta[k,othercategory]
    ccc[othercategory] <- 1
    cc <- ccc - ccc[K]
    b <- c(b, -ccc[K])
    A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  }
  
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  h <- rcdd::makeH(constr$constr, constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  polygon_df <- rbind(polygon_df, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], cat = k))
}

gpointsshaded <- gpoints +  geom_polygon(data = polygon_df %>% filter(cat == 1), alpha = 0.5, aes(fill = factor(cat)))
gpointsshaded <- gpointsshaded + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpointsshaded

# ggsave(filename = "stepbystep.feasibleset.3.pdf", plot = gpointsshaded, width = 5, height = 5)

gpointsshaded <- gpoints +  geom_polygon(data = polygon_df %>% filter(cat <= 2), alpha = 0.5, aes(fill = factor(cat))) 
gpointsshaded <- gpointsshaded + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpointsshaded

# ggsave(filename = "stepbystep.feasibleset.4.pdf", width = 5, height = 5, plot = gpointsshaded)

gpointsshaded <- gpoints +  geom_polygon(data = polygon_df %>% filter(cat <= 3), alpha = 0.5, aes(fill = factor(cat)))
gpointsshaded <- gpointsshaded + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpointsshaded

# ggsave(filename = "stepbystep.feasibleset.5.pdf", width = 5, height = 5, plot = gpointsshaded)

## feasible set 

A <- Asimplex
b <- bsimplex
for (k in 1:K){
  for (othercategory in setdiff(1:K, k)){
    # then we add the constraint theta_k < b
    ccc <- rep(0, K)
    ccc[k] <- -eta[k,othercategory]
    ccc[othercategory] <- 1
    cc <- ccc - ccc[K]
    b <- c(b, -ccc[K])
    A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  }
}

constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
h <- rcdd::makeH(constr$constr, constr$rhs)
## try to find V representation (for Vendetta)
v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
  stop("Failed to enumerate vertices. Is the polytope unbounded?")
}
vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
# then add last coordinate, so that entries sum to one again
vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
average_ <- colMeans(vertices_cart)
o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
vertices_cart <- vertices_cart[o_,]

gpointsshaded <- gpoints + geom_polygon(data = data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
gpointsshaded <- gpointsshaded + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpointsshaded

# ggsave(filename = "stepbystep.feasibleset.6.pdf", width = 5, height = 5, plot = gpointsshaded)

## if we take any point in the feasible set ...
v1 <- graphsettings$v_cartesian[[1]]
v2 <- graphsettings$v_cartesian[[2]]
v3 <- graphsettings$v_cartesian[[3]]

pt_xy <- colMeans(vertices_cart)
gpoints_split <- gpoints + geom_point(aes(x = pt_xy[1], y = pt_xy[2])) + geom_polygon(data = data.frame(x = vertices_cart[,1], y= vertices_cart[,2]), alpha = 0.5)
gpoints_split <- gpoints_split + geom_polygon(data = data.frame(x = c(pt_xy[1], v2[1], v3[1]), y = c(pt_xy[2], v2[2], v3[2])),
                      colour = "black", alpha = 0.5)
gpoints_split <- gpoints_split + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v3[1]), y = c(pt_xy[2], v1[2], v3[2])),
                      colour = "black", alpha = 0.5)
gpoints_split <- gpoints_split + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v2[1]), y = c(pt_xy[2], v1[2], v2[2])),
                      colour = "black", alpha = 0.5)
gpoints_split <- gpoints_split + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpoints_split
# ggsave(filename = "stepbystep.feasibleset.7.pdf", width = 5, height = 5, plot = gpoints_split)

## now with a non feasible 
index <- which(!accept)[2]
eta <- etas[[index]]
pts_barcoord <- pts[[index]]

pts_cart <- lapply(pts_barcoord, function(l) t(apply(matrix(l, ncol = K), 1, function(v) barycentric2cartesian(v, v_cartesian))))

## draw some points uniformly on the simplex

df <- data.frame()
for (d in categories){
  df <- rbind(df, data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2], category = d))
}
df

gpoints <- g + geom_point(data=df, aes(x = x, y = y, fill = factor(category), col = factor(category)), size = 3, shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpoints

# ggsave(filename = "anothertrial.pdf", width = 5, height = 5, plot = gpoints)


polygon_df <- data.frame()
for (k in 1:K){
  A <- Asimplex
  b <- bsimplex
  
  for (othercategory in setdiff(1:K, k)){
    # then we add the constraint theta_k < b
    ccc <- rep(0, K)
    ccc[k] <- -eta[k,othercategory]
    ccc[othercategory] <- 1
    cc <- ccc - ccc[K]
    b <- c(b, -ccc[K])
    A <- rbind(A, matrix(cc[1:(K-1)], nrow = 1))
  }
  
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  h <- rcdd::makeH(constr$constr, constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cart <- t(apply(vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  average_ <- colMeans(vertices_cart)
  o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cart <- vertices_cart[o_,]
  polygon_df <- rbind(polygon_df, data.frame(x = vertices_cart[,1], y= vertices_cart[,2], cat = k))
}

gpointsshaded <- g + geom_polygon(data = polygon_df %>% filter(cat <= 3), alpha = 0.5, aes(fill = factor(cat))) + geom_point(data=df, aes(x = x, y = y, fill = factor(category), col = factor(category)), size = 3, shape = 21)
gpointsshaded <- gpointsshaded + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpointsshaded


# ggsave(filename = "infeasibleset.pdf", width = 5, height = 5, plot = gpointsshaded)
##


mean(accept) * 100
