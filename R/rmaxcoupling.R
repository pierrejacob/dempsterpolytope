## sample from maximum coupling of two distribution,
## one is uniform on Delta_k(theta_star1)
## the other is uniform on Delta_k(theta_star2)
## outputs the pair of samples, and an indicator that they are equal
#'@export
rmaxcoupling <- function(k, theta_star1, theta_star2){
  x <- montecarlodsm:::runif_piktheta_one_cpp(k, theta_star1)
  pdf1_x <- 1/theta_star1[k]
  pdf2_x <- montecarlodsm:::dunif_piktheta_cpp(x, k, theta_star2)
  if (runif(1) < (pdf2_x / pdf1_x)){
    return(list(pts = cbind(x,x), equal = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- montecarlodsm:::runif_piktheta_one_cpp(k, theta_star2)
      pdf2_y <- 1/theta_star2[k]
      pdf1_y <- montecarlodsm:::dunif_piktheta_cpp(y, k, theta_star1)
      reject <- (runif(1) < (pdf1_y/pdf2_y))
    }
    return(list(pts = cbind(x,y), equal = FALSE))
  }
}
