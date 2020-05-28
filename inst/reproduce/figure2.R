set.seed(1)
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
labelsize <- 5
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.65, y = -0.1, label = TeX("$\\theta_3/\\theta_2 = \\eta_{2\\rightarrow 3}$", output = "character")), 
                    aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.2, y = 0.6, label = TeX("$\\theta_1/\\theta_2 = \\eta_{2\\rightarrow 1}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.3, y = -0.1, label = TeX("$\\theta_2/\\theta_3 = \\eta_{3\\rightarrow 2}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[3], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.9, y = 0.4, label = TeX("$\\theta_1/\\theta_3 = \\eta_{3\\rightarrow 1}$", output = "character")), 
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
print(gcond)
ggsave(filename = "conditional2.pdf", plot = gcond, width = 5, height = 5)



