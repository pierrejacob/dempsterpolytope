rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_custom_theme()
set.seed(1)


counts <- c(8,4)
nsamples <- 10
exact_dirichlet_pts <- gtools::rdirichlet(nsamples, counts+c(2,2))
cvxpolytope_cartesian.df <- data.frame()
for (isample in 1:nsamples){
  cvx_Vrepresentation <- rbind(c(exact_dirichlet_pts[isample,], 0), c(0,0,1))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, isample = isample))
}
## 
g <- create_plot_triangle(graphsettings = graphsettings)
gdirichlet <- g + geom_path(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = isample), alpha = .25, lineend = 'round')
gdirichlet
##
ggsave(plot = gdirichlet, filename = "pp.priorsegments.pdf", width = 5, height = 5)

K <- 3
new_counts <- c(2,1,3)
new_niterations <- 1e3
new_results <- gibbs_sampler(new_niterations, new_counts)
subiter <- floor(seq(from = 1e2, to = new_niterations, length.out = 10))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  vert <- etas_vertices(new_results$etas[iter,,])
  cvx_cartesian <- t(apply(vert, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gpolytopes <- g + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 0.25, alpha = .25, fill = 'black', colour = 'black')
gpolytopes
ggsave(plot = gpolytopes, filename = "pp.polytopes.pdf", width = 5, height = 5)

gdirichletpolytopes <- gdirichlet + 
  geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), alpha = .25, size = 0.25, fill = 'black', colour = 'black')
gdirichletpolytopes

ggsave(plot = gdirichletpolytopes, filename = "pp.overlaid.pdf", width = 5, height = 5)

## intersection segment - polytope

## get constr representation of segment linking vertex 3 to point on edge1-2
segment_constr <- function(point){
  K_ <- 3
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K_-1, K_-1))
  b <- c(1, rep(0, K_-1))
  ccc <- c(-point[2]/point[1], 1, 0)
  cc <- ccc - ccc[K_]
  A <- rbind(A, matrix(cc[1:(K_-1)], nrow = 1))
  b <- c(b, -ccc[K_])
  ccc <- c(point[2]/point[1], -1, 0)
  cc <- ccc - ccc[K_]
  A <- rbind(A, matrix(cc[1:(K_-1)], nrow = 1))
  b <- c(b, -ccc[K_])
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  return(constr)
}

## for a convex polytope, draw segment and find intersection
intersect_segment <- function(eta){
  K_ <- dim(eta)[1]
  categories <- 1:K_
  # # the constraints are on the first K-1 coordinates
  # # the first ones say that the feasible set is within the simplex
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K_-1, K_-1))
  b <- c(1, rep(0, K_-1))
  # then the extra constraints come from etas
  for (d in categories){
    for (j in setdiff(categories, d)){
      if (is.finite(eta[d,j])){
        # cccc (wA wB wC ... )' = 0
        ccc <- rep(0, K_)
        ccc[d] <- -eta[d, j]
        ccc[j] <- 1
        cc <- ccc - ccc[K_]
        b <- c(b, -ccc[K_])
        A <- rbind(A, matrix(cc[1:(K_-1)], nrow = 1))
      } else {
        # if eta is infinite, no constraint
      }
    }
  }
  polytope_constr_ <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  
  segment_constr_ <- segment_constr(gtools::rdirichlet(1, counts+c(2,2)))
  constr <- rbind(polytope_constr_$constr, segment_constr_$constr)
  rhs <- c(polytope_constr_$rhs, segment_constr_$rhs)
  h <- rcdd::makeH(constr, rhs)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  if (dim(vertices_barcoord)[1] == 0){
    ## no intersection
    cvx_cartesian <- NULL
  } else {
    cvx_cartesian <- t(apply(vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
    average_ <- colMeans(cvx_cartesian)
    o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
    cvx_cartesian <- cvx_cartesian[o_,]
  }
  return(cvx_cartesian)
}

## loop over Gibbs output and intersect with some random segments
subiter <- 1e2:new_niterations

intersections <- lapply(subiter, function(i) intersect_segment(new_results$etas[i,,]))
mean(sapply(intersections, function(x) is.null(x)))
intersections.df <- lapply(1:length(intersections), function(i){
  if (!is.null(intersections[[i]])){
    data.frame(x = intersections[[i]][,1], y = intersections[[i]][,2], iter = i)
  } else {
    data.frame()
  }
}) %>% bind_rows()

(intersections.df %>% pull(iter) %>% unique)[10]
gintersection <- g + geom_polygon(data = intersections.df %>% filter(iter <= 200), aes(x = x, y = y, group = iter),
                                  size = 0.25, alpha = .25, fill = 'black', colour = 'black')
gintersection

gintersection <- g + geom_path(data = intersections.df %>% filter(iter <= 200), aes(x = x, y = y, group = iter),
                               alpha = .25, lineend = "round")


ggsave(plot = gintersection, filename = "pp.intersection.pdf", width = 5, height = 5)

