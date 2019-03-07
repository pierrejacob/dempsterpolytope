
#'@export
set_my_theme <- function(){
  library(ggplot2)
  library(ggthemes)
  library(gridExtra)  
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 20, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 20, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               legend.position = "bottom",
               plot.margin = unit(c(1,1,1,1), "cm"))
}

# Convert barycentric coordinates to cartesian coordinates
#'@export
barycentric2cartesian <- function(barycentric, v_cartesian){
  return(rowSums(sapply(1:length(v_cartesian), function(i) barycentric[i] * v_cartesian[[i]])))
}
# Convert back 
#'@export
cartesian2barycentric <- function(cartesian, matrixT, v_cartesian){
  w <- solve(matrixT, matrix(cartesian - v_cartesian[[length(v_cartesian)]], ncol = 1))[,1]
  return(c(w, 1 - sum(w)))
}

# get intersection of two lines, where "interslopei" is 
# a vector containing the intercept and the slope of line i
# (all in x,y coordinates)
#'@export
get_line_intersection <- function(interslope1, interslope2){
  xinter <- (interslope1[1] - interslope2[1])/(interslope2[2] - interslope1[2])
  yinter <- interslope1[1] + xinter * interslope1[2]
  return(c(xinter, yinter))
}

# function to plot triangle with ggplot
#'@export
ggplot_triangle <- function(v_cartesian, pts_cartesian, etas, addpolytope = FALSE, removelines = FALSE, cols = c("red", "green", "blue")){
  if (missing(pts_cartesian)){
    addpoints <- FALSE
  } else {
    addpoints <- TRUE
  }
  if (missing(etas)){
    addconstraints <- FALSE
  } else {
    addconstraints <- TRUE
  }
  K <- 3
  categories <- 1:K
  matrixT <- matrix(0, nrow = K-1, ncol = K-1)
  for (k in 1:(K-1)){
    kth_components <- sapply(v_cartesian, function(x) x[k])
    matrixT[k,] <- kth_components[-K] - kth_components[K]
  }
  v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
  triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
  g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
  g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
  g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
  g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
  g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
  if (addpoints){
    for (d in categories){
      g <- g + geom_point(data=data.frame(x = pts_cartesian[[d]][,1], y = pts_cartesian[[d]][,2]), colour = cols[d])
    }
    g <- g + scale_color_discrete("category: ")
  }
  if (addconstraints){
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
    for (d in categories){
      # set indices for two other components
      j1 <- setdiff(categories, d)[1]
      j2 <- setdiff(categories, d)[2]
      interslope_j1 <- barconstraint2cartconstraint(d, j1, 1/etas[d, j1], matrixT, v_cartesian)
      interslope_j2 <- barconstraint2cartconstraint(d, j2, 1/etas[d, j2], matrixT, v_cartesian)
      if (!(removelines)){
        g <- g + geom_abline(intercept = interslope_j1[1], slope = interslope_j1[2], colour = cols[d], linetype = 2)
        g <- g + geom_abline(intercept = interslope_j2[1], slope = interslope_j2[2], colour = cols[d], linetype = 2)
      }
      intersection_12 <- get_line_intersection(interslope_j1, interslope_j2)
      g <- g + geom_polygon(data = data.frame(x = c(v_cartesian[[j1]][1], v_cartesian[[j2]][1], intersection_12[1]),
                                              y = c(v_cartesian[[j1]][2], v_cartesian[[j2]][2], intersection_12[2])), alpha = 0.2, fill = cols[d])
    }
    if (addpolytope){
      etascvxp <- etas2cvxpolytope(etas)
      ## convert coordinates to cartesian
      vertices_cart <- t(apply(etascvxp$vertices_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
      # order vertices according to angles
      average_ <- colMeans(vertices_cart)
      o_ <- order(apply(sweep(vertices_cart, 2, average_, "-"), 1, function(v) atan2(v[2], v[1])))
      vertices_cart <- vertices_cart[o_,]
      g <- g + geom_polygon(data=data.frame(x = vertices_cart[,1], y= vertices_cart[,2]))
    }
  }
  return(g)
}

