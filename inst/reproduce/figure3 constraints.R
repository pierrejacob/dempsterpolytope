set.seed(1)
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
attach(graphsettings)
K <- 3
categories <- 1:K

counts <- c(2,3,1)
niterations <- 100

gibbs_results <- gibbs_sampler(niterations, counts)

# function to plot triangle with ggplot
g <- create_plot_triangle(graphsettings)

iter <- 70
dim(gibbs_results$Us[[1]])

gconstraints <- g
eta <- gibbs_results$etas[iter,,]
for (d in categories){
  for (j in setdiff(categories, d)) gconstraints <- add_ratioconstraint(graphsettings, gconstraints, eta[d,j], d, j, graphsettings$contcols[d])
}
labelsize <- 5
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.65, y = -0.1, label = TeX("$\\theta_3/\\theta_2 = \\eta_{2\\rightarrow 3}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.2, y = 0.6, label = TeX("$\\theta_1/\\theta_2 = \\eta_{2\\rightarrow 1}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[2], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.3, y = -0.1, label = TeX("$\\theta_2/\\theta_3 = \\eta_{3\\rightarrow 2}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[3], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.9, y = 0.4, label = TeX("$\\theta_1/\\theta_3 = \\eta_{3\\rightarrow 1}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[3], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.1, y = 0.2, label = TeX("$\\theta_2/\\theta_1 = \\eta_{1\\rightarrow 2}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[1], alpha = 0.75, size = labelsize)
gconstraints <- gconstraints + geom_label(data = data.frame(x = 0.92, y = 0.18, label = TeX("$\\theta_3/\\theta_1 = \\eta_{1\\rightarrow 3}$", output = "character")), 
                                          aes(x = x, y = y, label = label), parse = TRUE, col = contcols[1], alpha = 0.75, size = labelsize)
gconstraints
ggsave(filename = "etaconstraints.pdf", plot = gconstraints, width = 5, height = 5)
# 

library(igraph)
library(ggraph)

# counts <- c(2,3,1)
# niterations <- 100
# gibbs_results <- gibbs_sampler(niterations, counts)
# eta <- gibbs_results$etas[niterations,,]

graph_ <- igraph::graph_from_adjacency_matrix(adjmatrix = matrix(1, nrow = 3, ncol = 3), mode = "directed", diag = FALSE)

V(graph_)$name = c(1,2,3)
E(graph_)$name <- paste(E(graph_))
labelsize <- 5

g <- ggraph(graph_, x = sapply(v_cartesian, function(x) x[1]), y = sapply(v_cartesian, function(x) x[2])) 
g <- g + geom_edge_fan(arrow = arrow(), end_cap = circle(7, 'mm')) +  theme_graph(fg_text_colour = 'white')
g <- g + geom_node_label(aes(label = name), size = 10, fill = "white")
g <- g + annotate(geom = "label", x = 0.15, y = 0.5, label = TeX("$\\log(\\eta_{1\\rightarrow 2})$"), size = labelsize, col = contcols[1])
g <- g + annotate(geom = "label", x = 0.35, y = 0.4, label = TeX("$\\log(\\eta_{2\\rightarrow 1})$"), size = labelsize, col = contcols[2])
g <- g + annotate(geom = "label", x = 0.85, y = 0.5, label = TeX("$\\log(\\eta_{3\\rightarrow 1})$"), size = labelsize, col = contcols[3])
g <- g + annotate(geom = "label", x = 0.65, y = 0.4, label = TeX("$\\log(\\eta_{1\\rightarrow 3})$"), size = labelsize, col = contcols[1])
g <- g + annotate(geom = "label", x = 0.5, y = +0.1, label = TeX("$\\log(\\eta_{3\\rightarrow 2})$"), size = labelsize, col = contcols[3])
g <- g + annotate(geom = "label", x = 0.5, y = -0.1, label = TeX("$\\log(\\eta_{2\\rightarrow 3})$"), size = labelsize, col = contcols[2])
g 
# 
ggsave(filename = "fullyconnectedgraph.pdf", plot = g, width = 5, height = 5)

