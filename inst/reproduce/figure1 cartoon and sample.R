## This script generates the two plots of Figure 1 of the article.
rm(list = ls())
library(dempsterpolytope)
library(latex2exp)
graphsettings <- set_custom_theme()
set.seed(1)

# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(2,3,1)

## triangle with equal sides
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
pt_bar <- c(0.3, 0.25, 0.45)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)
# Now with ggplot2
g <- create_plot_triangle(graphsettings)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 1, fill = graphsettings$cols[1], alpha = 0.8)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 2, fill = graphsettings$cols[2], alpha = 0.8)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 3, fill = graphsettings$cols[3], alpha = 0.8)
meanpi1 <- c(pt_xy[1]/3 + v2[1]/3 + v3[1]/3, pt_xy[2]/3 + v2[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(Delta[1](theta)), parse = FALSE, size = 5) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(Delta[2](theta)), parse = FALSE, size = 5) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(Delta[3](theta)), parse = FALSE, size = 5) 
g <- g + annotate(geom = 'text', size = 5, x = .56, y = .28, label = TeX("$\\theta$", output = "character"), parse = TRUE)
g
ggsave(filename = "subsimplices.pdf", plot = g, width = 5, height = 5)


## create subsimplex as convex polytope
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

ntrials <- 2e2
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
index <- which(accept)[2]
eta <- etas[[index]]
pts_barcoord <- pts[[index]]

## draw some points uniformly on the simplex
g <- create_plot_triangle(graphsettings)
# ## points are each related to one of the x
gpoints <- g
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[1]], colour = graphsettings$contcols[1], fill = graphsettings$cols[1])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[2]], colour = graphsettings$contcols[2], fill = graphsettings$cols[2])
gpoints <- add_plot_points(graphsettings, gpoints, pts_barcoord[[3]], colour = graphsettings$contcols[3], fill = graphsettings$cols[3])

attach(graphsettings)
eta_cvx <- etas2cvxpolytope(eta)
gpoints_polytope <- add_plot_polytope(graphsettings, gpoints, eta_cvx)
gpoints_polytope
##
ggsave(filename = "onesample.pdf", plot = gpoints_polytope, width = 5, height = 5)

