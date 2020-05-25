library(igraph)
library(ggraph)
library(dempsterpolytope)
library(latex2exp)
graphsettings <- set_custom_theme()
set.seed(1)
# counts <- c(2,3,1)
# niterations <- 100
# gibbs_results <- gibbs_sampler(niterations, counts)
# eta <- gibbs_results$etas[niterations,,]

graph_ <- igraph::graph_from_adjacency_matrix(adjmatrix = matrix(1, nrow = 3, ncol = 3), mode = "directed", diag = FALSE)

V(graph_)$name = c(1,2,3)
E(graph_)$name <- paste(E(graph_))
g <- ggraph(graph_, layout = "kk") + geom_edge_fan(arrow = arrow(), end_cap = circle(7, 'mm')) +  theme_graph(fg_text_colour = 'white')
g <- g + geom_node_label(aes(label = name), size = 10)
g <- g + annotate(geom = "text", x = 0.5, y = 0.6, label = TeX("$\\log(\\eta_{1\\rightarrow 2})$"), size = 8)
g <- g + annotate(geom = "text", x = 0.1, y = 0.20, label = TeX("$\\log(\\eta_{2\\rightarrow 1})$"), size = 8)
g

# ggsave(filename = "fullyconnectedgraph.pdf", plot = g, width = 7, height = 7)
