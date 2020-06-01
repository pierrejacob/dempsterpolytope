# devtools::install_github("pierrejacob/dempsterpolytope")
rm(list = ls())
library(dempsterpolytope)
library(doParallel)
library(doRNG)
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
registerDoParallel(cores = detectCores()-2)
graphsettings <- set_my_theme()
set.seed(1)

names(graphsettings)

gtriangle <- dempsterpolytope::ggplot_triangle(graphsettings$v_cartesian)

etas <- matrix(Inf, 3, 3)
etas[1,2] <- 5
etas[2,1] <- 1/3
cvx <- etas2cvxpolytope(etas)
cvx
dempsterpolytope:::check_cst_graph(etas)

## plot convex polytope as polygon 
cvx2cart <- function(cvx, graphsettings){
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  return(data.frame(x = cvx_cartesian[,1], y = cvx_cartesian[,2]))
}

# gtriangle + geom_polygon_pattern(data=cvx2cart(cvx, graphsettings), aes(x = x, y = y),
#                                  pattern_spacing = .05, pattern_density = .5, 
#                                  pattern_angle = 90, pattern_fill = "white", fill = 'black')

g <- gtriangle + geom_polygon_pattern(data=cvx2cart(cvx, graphsettings), aes(x = x, y = y),
                                 pattern_spacing = .05, pattern_density = .5, 
                                 pattern_angle = 90, pattern_fill = "white", fill = 'black')

## the constraint encodes theta_d/theta_j = eta
## the output is an intercept and a slope in cartesian coordinate
barconstraint2cartconstraint <- function(d, j, eta, graphsettings){
  # ccc * (wA wB wC) = 0 with:
  ccc <- rep(0, 3)
  ccc[d] <- 1
  ccc[j] <- - eta
  # which is equivalent to ftilde * (wA wB) = gtilde with
  ftilde <- c(ccc[1] - ccc[3], ccc[2] - ccc[3])
  gtilde <- -ccc[3]
  # we can generically express that as a constraint on x,y through
  f <- solve(t(graphsettings$matrixT), ftilde)
  g <- gtilde + sum(f * graphsettings$v_cartesian[[3]])
  # f1 x + f2 y = g is equivalent to a = g/f2, b = - f1/f2
  return(c(g/f[2], -f[1]/f[2]))
}

segment_thetad_over_thetaj_equal_eta <- function(d, j, eta, graphsettings){
  K_ <- 3
  A <- matrix(rep(1, K_-1), ncol = K_-1)
  A <- rbind(A, diag(-1, K_-1, K_-1))
  b <- c(1, rep(0, K_-1))
  # ccc * (wA wB wC) = 0 with:
  ccc <- rep(0, 3)
  ccc[d] <- 1
  ccc[j] <- - eta
  # which is equivalent to ftilde * (wA wB) = gtilde with
  ftilde <- c(ccc[1] - ccc[3], ccc[2] - ccc[3])
  gtilde <- -ccc[3]
  A <- rbind(A, ftilde)
  b <- c(b, gtilde)
  A <- rbind(A, -ftilde)
  b <- c(b, -gtilde)
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  ## make H representation
  h <- rcdd::makeH(constr$constr, constr$rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cartesian <- t(apply(vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
  average_ <- colMeans(vertices_cartesian)
  o_ <- order(apply(sweep(vertices_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  vertices_cartesian <- vertices_cartesian[o_,]
  return(list(vertices_barcoord = vertices_barcoord, constr = constr, vertices_cartesian = vertices_cartesian))
}

# # # segment_thetad_over_thetaj_equal_eta <- function(d, j, eta, graphsettings){
# # d <- 1
# # j <- 2 
# # eta <- etas[2,1]
# # eta_tilde <- eta 
# # K_ <- 3
# # A <- matrix(rep(1, K_-1), ncol = K_-1)
# # A <- rbind(A, diag(-1, K_-1, K_-1))
# # b <- c(1, rep(0, K_-1))
# # # ccc * (wA wB wC) = 0 with:
# # ccc <- rep(0, 3)
# # ccc[d] <- 1
# # ccc[j] <- - eta
# # # which is equivalent to ftilde * (wA wB) = gtilde with
# # ftilde <- c(ccc[1] - ccc[3], ccc[2] - ccc[3])
# # gtilde <- -ccc[3]
# # A <- rbind(A, ftilde)
# # b <- c(b, gtilde)
# # A <- rbind(A, -ftilde)
# # b <- c(b, -gtilde)
# # constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
# # ## make H representation
# # h <- rcdd::makeH(constr$constr, constr$rhs)
# # ## try to find V representation (for Vendetta)
# # v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
# # if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
# #   stop("Failed to enumerate vertices. Is the polytope unbounded?")
# # }
# # vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
# # # then add last coordinate, so that entries sum to one again
# # vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
# # vertices_barcoord
# 
# 
# 
# vertices_cartesian <- t(apply(vertices_barcoord, 1, function(row) barycentric2cartesian(row, graphsettings$v_cartesian)))
# average_ <- colMeans(vertices_cartesian)
# o_ <- order(apply(sweep(vertices_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
# vertices_cartesian <- vertices_cartesian[o_,]
# # return(arrows_cartesian))
# # }



abline12 <- barconstraint2cartconstraint(1, 2, etas[2,1], graphsettings)
abline21 <- barconstraint2cartconstraint(2, 1, etas[1,2], graphsettings)
gtriangle + geom_abline(intercept = abline12[1], slope = abline12[2], linetype = 3) + 
  geom_abline(intercept = abline21[1], slope = abline21[2])

seg1 <- segment_thetad_over_thetaj_equal_eta(2, 1, etas[1,2], graphsettings)$vertices_cartesian
gtriangle + geom_path(data = data.frame(seg1), aes(x = X1, y = X2))




