library(dempsterpolytope)
library(doParallel)
library(doRNG)
library(latex2exp)
registerDoParallel(cores = detectCores()-2)
set_my_theme()
set.seed(1)
rm(list = ls())
## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
contcols <- c("darkred", "goldenrod3", "darkblue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
set.seed(4)
##
K <- 3
categories <- 1:K
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
gtriangle <- g
##

K <- 2
counts <- c(3,2)
niterations <- 1e3
results <- gibbs_sampler(niterations, counts)

##
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 20))
offset <- truncnorm::rtruncnorm(n = niterations, a = -.02, b = .02, sd = .01)
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(results$etas[iter,,])
  cvx_Vrepresentation <- cbind(cvx$vertices_barcoord, 0)
  cvx_Vrepresentation <- rbind(cvx_Vrepresentation,
                       cbind(cvx$vertices_barcoord, offset[iter]))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
gintervals <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                      alpha = .4)
# ggsave(plot = gintervals, filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/minext1.pdf",
#        width = 5, height = 5)

cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  etas <- matrix(1, nrow = 3, ncol = 3)
  etas[1:2,1:2] <- results$etas[iter,,]
  etas[1,3] <- etas[2,3] <- Inf
  etas[3,1] <- etas[3,2] <- Inf
  cvx <- etas2cvxpolytope(etas)
  cvx_Vrepresentation <- cvx$vertices_barcoord
  cvx_Vrepresentation[1:2,3] <- offset[iter]
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}

gintervalext <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 2,
                      alpha = .4)
# ggsave(plot = gintervalext, filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/minext2.pdf",
#        width = 5, height = 5)

##
subiter <- floor(seq(from = 1e2, to = niterations, length.out = 30))
exact_dirichlet_pts <- gtools::rdirichlet(niterations, counts+c(2,2))
cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx_Vrepresentation <- rbind(c(exact_dirichlet_pts[iter,], 0), c(0,0,1))
  cvx_cartesian <- t(apply(cvx_Vrepresentation, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}

gdirichlet <- gtriangle + geom_path(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 1,
                                         alpha = .4)
ggsave(plot = gdirichlet, filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/minext3.png",
       width = 5, height = 5)

## 
K <- 3
new_counts <- c(0,0,3)
new_niterations <- 1e3
new_results <- gibbs_sampler(new_niterations, new_counts)

cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  cvx <- etas2cvxpolytope(new_results$etas[iter,,])
  cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  average_ <- colMeans(cvx_cartesian)
  o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
  cvx_cartesian <- cvx_cartesian[o_,]
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(cvx_cartesian, iter = iter))
}
g3 <- gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                                       alpha = .2)
ggsave(plot = g3, filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/minext4.png",
       width = 5, height = 5)


## intersection between polytope and segment
## polytope should be output of etas2cvxpolytope
## and segment should be a matrix with two rows, each row giving vertex in barycentric coordinates
iter <- subiter[3]
cvx <- etas2cvxpolytope(new_results$etas[iter,,])
cvx_cartesian <- t(apply(cvx$vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
average_ <- colMeans(cvx_cartesian)
o_ <- order(apply(sweep(cvx_cartesian, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
cvx_cartesian <- cvx_cartesian[o_,]
cvxpolytope_cartesian.df <- data.frame(cvx_cartesian, iter = iter)

segment <- rbind(c(exact_dirichlet_pts[iter,], 0), c(0,0,1))
segment <- t(apply(segment, 1, function(row) barycentric2cartesian(row, v_cartesian)))
segment.df <- data.frame(segment, iter = iter)

gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                         alpha = .4) + geom_path(data = segment.df, aes(x = X1, y = X2, group = iter))

intersect_polytope_ <- function(polytope, pt_12){
  lineconstr <- c(1,-pt_12[1]/pt_12[2])
  lineconstr <- rbind(lineconstr, c(-1,pt_12[1]/pt_12[2]))
  constr <- rbind(polytope$constr$constr, lineconstr)
  rhs <- c(polytope$constr$rhs, c(0,0))
  dir <- c(polytope$constr$dir, c('<=', '<='))
  h <- rcdd::makeH(constr, rhs)
  ## try to find V representation (for Vendetta)
  v <- rcdd::q2d(rcdd::scdd(rcdd::d2q(h))$output)
  if (any(v[, 1] != "0") || any(v[, 2] != "1")) {
    stop("Failed to enumerate vertices. Is the polytope unbounded?")
  }
  vertices_barcoord <- v[, -c(1, 2), drop = FALSE]
  # then add last coordinate, so that entries sum to one again
  vertices_barcoord <- cbind(vertices_barcoord, 1- apply(vertices_barcoord, 1, sum))
  vertices_cart <- t(apply(vertices_barcoord, 1, function(row) barycentric2cartesian(row, v_cartesian)))
  return(list(vertices_barcoord = vertices_barcoord, vertices_cart = vertices_cart))
}

intersect_cvx_ <- intersect_polytope_(etas2cvxpolytope(new_results$etas[iter,,]), exact_dirichlet_pts[iter,])


gtriangle + geom_polygon(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter),
                         alpha = .4) + geom_path(data = segment.df, aes(x = X1, y = X2, group = iter)) +
  geom_path(data=data.frame(intersect_cvx_$vertices_cart), aes(x = X1, y = X2), colour = 'orange')


cvxpolytope_cartesian.df <- data.frame()
for (iter in subiter){
  intersect_cvx_ <- intersect_polytope_(etas2cvxpolytope(new_results$etas[iter,,]), exact_dirichlet_pts[iter,])
  cvxpolytope_cartesian.df <- rbind(cvxpolytope_cartesian.df, data.frame(intersect_cvx_$vertices_cart, iter = iter))
}

gintersect <- gtriangle + geom_path(data = cvxpolytope_cartesian.df, aes(x = X1, y = X2, group = iter), size = 1,
                                    alpha = .4)

ggsave(plot = gintersect, filename = "~/Dropbox/Dempster-Shafer Analysis/GibbsSamplerMultinomial/minext5.png",
       width = 5, height = 5)

