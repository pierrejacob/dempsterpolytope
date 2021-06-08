set.seed(3)
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
attach(graphsettings)
##
K <- 3
categories <- 1:K
counts <- c(2,3,1)
## number of MCMC iterations
niterations <- 100
gibbs_results <- gibbs_sampler(niterations, counts)

# function to plot triangle with ggplot
g <- create_plot_triangle(graphsettings)

iter <- 70
dim(gibbs_results$Us[[1]])

gpoints <- add_plot_points(graphsettings, g = g, barypoints = gibbs_results$Us[[1]][iter,,], colour = graphsettings$contcols[1], fill = graphsettings$cols[1])
gpoints <- add_plot_points(graphsettings, g = gpoints, barypoints = gibbs_results$Us[[2]][iter,,], colour = graphsettings$contcols[2], fill = graphsettings$cols[2])
gpoints <- add_plot_points(graphsettings, g = gpoints, barypoints = gibbs_results$Us[[3]][iter,,], colour = graphsettings$contcols[3], fill = graphsettings$cols[3])
eta_vertices_ <- etas_vertices(gibbs_results$etas[iter,,])
gpoints <- add_plot_polytope(graphsettings, gpoints, eta_vertices_)
gpoints

ggsave(filename = "conditional1.pdf", plot = gpoints, width = 5, height = 5)

## now remove points of category 1
gpoints <- add_plot_points(graphsettings, g = g, barypoints = gibbs_results$Us[[2]][iter,,], colour = graphsettings$contcols[2], fill = graphsettings$cols[2])
gpoints <- add_plot_points(graphsettings, g = gpoints, barypoints = gibbs_results$Us[[3]][iter,,], colour = graphsettings$contcols[3], fill = graphsettings$cols[3])
gpoints


eta <- gibbs_results$etas[iter,,]
gconstraints <- gpoints
gconstraints <- add_ratioconstraint(graphsettings, gconstraints, eta[2,1], 2, 1, graphsettings$contcols[2])
gconstraints <- add_ratioconstraint(graphsettings, gconstraints, eta[2,3], 2, 3, graphsettings$contcols[2])
gconstraints <- add_ratioconstraint(graphsettings, gconstraints, eta[3,1], 3, 1, graphsettings$contcols[3])
gconstraints <- add_ratioconstraint(graphsettings, gconstraints, eta[3,2], 3, 2, graphsettings$contcols[3])
labelsize <- 7
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.6, y = -0.1, label = TeX("$\\theta_3/\\theta_2 = \\eta_{2\\rightarrow 3}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.15, y = 0.4, label = TeX("$\\theta_1/\\theta_2 = \\eta_{2\\rightarrow 1}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.15, y = -0.1, label = TeX("$\\theta_2/\\theta_3 = \\eta_{3\\rightarrow 2}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[3], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.85, y = 0.35, label = TeX("$\\theta_1/\\theta_3 = \\eta_{3\\rightarrow 1}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[3], alpha = 0.75, size = labelsize)

print(gconstraints)
# ggsave(filename = "conditional2.pdf", plot = gconstraints, width = 5, height = 5)

## find new conditional distribution for points in category 1
k <- 1
graph_ <- igraph::graph_from_adjacency_matrix(log(eta), mode = "directed", weighted = TRUE, diag = FALSE)
theta_star <- rep(0, K)
# minimum value among paths from k to ell
notk <- setdiff(1:K, k)
minimum_values <- rep(1, K)
# compute shortest paths in the graph
minimum_values[notk] <- igraph::distances(graph_, v = notk, to = k, mode = "out")[,1]
# compute solution based on values of shortest paths 
theta_star <- exp(-minimum_values)
theta_star[k] <- 1
theta_star <- theta_star / sum(theta_star)
gcond <- add_plot_subsimplex(graphsettings, gconstraints, theta_star, 1, fill = graphsettings$cols[1], alpha = 0.6)
gcond
theta_star
pts_k <- dempsterpolytope:::runif_piktheta_cpp(counts[k], k, theta_star)
pts_k
cartipoints <- t(apply(pts_k$pts, 1, function(v) barycentric2cartesian(v, graphsettings$v_cartesian)))
gcond <- gcond + geom_point(data=data.frame(x = cartipoints[,1], y = cartipoints[,2]), aes(x = x, y = y), colour = graphsettings$contcols[1], 
                   fill = graphsettings$cols[1], size = 3, shape = 15)

print(gcond)

ggsave(filename = "conditional2.pdf", plot = gcond, width = 5, height = 5)

### add feasible set
pts_k
g <- create_plot_triangle(graphsettings)
g <- add_plot_points(graphsettings, g = g, barypoints = pts_k$pts, colour = graphsettings$contcols[1], fill = graphsettings$cols[1])
g <- add_plot_points(graphsettings, g = g, barypoints = gibbs_results$Us[[2]][iter,,], colour = graphsettings$contcols[2], fill = graphsettings$cols[2])
g <- add_plot_points(graphsettings, g = g, barypoints = gibbs_results$Us[[3]][iter,,], colour = graphsettings$contcols[3], fill = graphsettings$cols[3])
# create etas
Us <- list()
Us[[1]] <- pts_k$pts
Us[[2]] <- gibbs_results$Us[[2]][iter,,]
Us[[3]] <- gibbs_results$Us[[3]][iter,,]
etas <- diag(1, 3, 3)
for (k in 1:K){
  notk <- setdiff(1:K, k)
  a_k <- matrix(Us[[k]], ncol = 3)
  for (ell in notk){
    etas[k,ell] <- min(a_k[,ell]/a_k[,k])
  }
}

eta_vert_ <- etas_vertices(etas)
gcond <- add_plot_polytope(graphsettings, gcond, eta_vert_)
gcond
ggsave(filename = "conditional3.pdf", plot = gcond, width = 5, height = 5)

