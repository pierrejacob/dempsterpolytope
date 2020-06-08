# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

## sample from maximum coupling of two distribution,
## one is uniform on Delta_k(theta_star1)
## the other is uniform on Delta_k(theta_star2)
## outputs the pair of samples, and an indicator that they are equal
#'@export
rmaxcoupling <- function(k, theta_star1, theta_star2){
  x <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star1)
  pdf1_x <- 1/theta_star1[k]
  pdf2_x <- dempsterpolytope:::dunif_piktheta_cpp(x, k, theta_star2)
  if (runif(1) < (pdf2_x / pdf1_x)){
    return(list(pts = cbind(x,x), equal = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- dempsterpolytope:::runif_piktheta_one_cpp(k, theta_star2)
      pdf2_y <- 1/theta_star2[k]
      pdf1_y <- dempsterpolytope:::dunif_piktheta_cpp(y, k, theta_star1)
      reject <- (runif(1) < (pdf1_y/pdf2_y))
    }
    return(list(pts = cbind(x,y), equal = FALSE))
  }
}
