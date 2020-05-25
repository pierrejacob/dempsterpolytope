library(dempsterpolytope)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
set_custom_theme()
set.seed(1)
rm(list = ls())

v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
##
K <- 3
categories <- 1:K
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}

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


counts <- c(9,8,3)
niterations <- 100

gibbs_results <- gibbs_sampler(niterations, counts)

# function to plot triangle with ggplot
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)


iter <- 70
dim(gibbs_results$Us[[1]])

pts_cart <- list()
for (d in categories){
  pts_cart[[d]] <- t(apply(gibbs_results$Us[[d]][iter,,], 1, function(v) barycentric2cartesian(v, v_cartesian)))
}

df <- data.frame()
for (d in categories){
  df <- rbind(df, data.frame(x = pts_cart[[d]][,1], y = pts_cart[[d]][,2], category = d))
}
df

gpoints <- g + geom_point(data=df, aes(x = x, y = y, fill = factor(category), col = factor(category)), shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpoints
# ggsave(filename = "conditional.1.pdf", plot = gpoints, width = 5, height = 5)
## now remove points of category 1
gpoints <- g + geom_point(data=df %>% filter(category != 1), aes(x = x, y = y, fill = factor(category), col = factor(category)), shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols[2:3]) + scale_color_manual("", values = contcols[2:3]) + theme(legend.position = "none")
gpoints
# ggsave(filename = "conditional.2.pdf", plot = gpoints, width = 5, height = 5)

eta <- gibbs_results$etas[iter,,]
for (d in setdiff(categories, 1)){
  # set indices for two other components
  j1 <- setdiff(categories, d)[1]
  j2 <- setdiff(categories, d)[2]
  interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/eta[d, j1], matrixT, v_cartesian)
  interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/eta[d, j2], matrixT, v_cartesian)
  g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = contcols[d], linetype = 2)
  g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = contcols[d], linetype = 2)
  intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
}
gpoints <- g + geom_point(data=df%>% filter(category != 1), aes(x = x, y = y, fill = factor(category), col = factor(category)), shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols[2:3]) + scale_color_manual("", values = contcols[2:3]) + theme(legend.position = "none")
gpoints
# ggsave(filename = "conditional.3.pdf", plot = gpoints, width = 5, height = 5)

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
theta_star
theta_star_cart <- barycentric2cartesian(theta_star, v_cartesian)
gcond <- gpoints + geom_polygon(data = data.frame(x = c(theta_star_cart[1], v2[1], v3[1]), y = c(theta_star_cart[2], v2[2], v3[2])), 
                       fill = "red", alpha = 0.5)
gcond
# ggsave(filename = "conditional.4.pdf", plot = gcond, width = 5, height = 5)



