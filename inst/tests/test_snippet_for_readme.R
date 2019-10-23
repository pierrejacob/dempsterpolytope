library(dempsterpolytope)
set.seed(1)

# count data
counts <- c(12, 14, 7)
# number of MCMC iterations
niterations <- 1e2
# run Gibbs sampler
gibbs_results <- gibbs_sampler(niterations, counts)

# obtain K x K matrix 
eta <- gibbs_results$etas_chain[niterations,,]

# convert it to H-representation and V-representation of a convex polytope
eta_converted <- etas2cvxpolytope(eta)
# H-representation
eta_converted$constr
# V-representation
eta_converted$vertices_barcoord

# next we can view the polytope in the K-simplex as a triangle, since here K = 3

v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
## this loads ggplot2
set_my_theme()
g <- ggplot_triangle(v_cartesian, etas = eta, addpolytope = T, cols = cols)
g
