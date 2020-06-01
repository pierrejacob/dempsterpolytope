## 
rm(list = ls())
library(dempsterpolytope)
graphsettings <- set_custom_theme()

x <- seq(-1,1,0.01) 
y <- seq(-1,1,0.01)
f <- function(x,y){ z <- -x - y + 1 }
z <- outer(x,y,f)
z <- ifelse(z<0,NA,z)
# persp_res <- persp(x, y, z, xlim = c(.4,1), ylim = c(0.4, 0.445), theta = 30, phi = 30, expand = 0.5,
#                    r = 1, col = "white", border = NA, box = FALSE, axes = 0)

v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))

# v1 <- c(0,0,0)
# v2 <- c(1,0,0)
# v3 <- c(1/2,sin(pi/3),0)
# v4 <- c(1/2,sin(pi/3)/2,1)

v1 <- c(sqrt(8/9),0,-1/3)
v2 <- c(-sqrt(2/9),sqrt(2/3),-1/3)
v3 <- c(-sqrt(2/9),-sqrt(2/3),-1/3)
v4 <- c(0,0,1)
z <- matrix(rep(-1/3, 4), ncol = 2)
z[1,] <- 1
persp_res <- persp(x = c(-sqrt(2/9), sqrt(8/9)), y = c(-sqrt(2/3), sqrt(2/3)), z = z, 
                   r = 2, d = 2,
                   xlim = c(-.5,.6), zlim = c(-.45,0.95),
                   theta = 200, phi = 1, expand = 1, col = "white", border = NA, box = FALSE, axes = 0)

add_segment <- function(v, vprime, persp_res){
  seg_ <- cbind(v, vprime)
  lines(trans3d(seg_[1,], seg_[2,], seg_[3,], pmat = persp_res), col = rgb(0,0,0))
}



add_segment(v1, v2, persp_res)
add_segment(v1, v3, persp_res)
add_segment(v2, v3, persp_res)
add_segment(v1, v4, persp_res)
add_segment(v2, v4, persp_res)
add_segment(v3, v4, persp_res)

## 
u_grid <- rev(seq(from = 0, to = 1, by = 0.02))
for (iu in 2:length(u_grid)){
  uprev <- u_grid[iu-1]
  unext <- u_grid[iu]
  polyg1 <- v1 * uprev + v2 * (1-uprev)
  polyg2 <- v1 * unext + v2 * (1-unext)
  polyg3 <- v3 * unext + v4 * (1-unext)
  polyg4 <- v3 * uprev + v4 * (1-uprev)
  polyg <- rbind(polyg1, polyg2, polyg3, polyg4)
  # add_segment(v_14, v_23, persp_res)
  polygon(trans3d(polyg[,1], polyg[,2], polyg[,3], pmat = persp_res), col = 'white', border = 'black')
}

text(trans3d(v1[1]+0.1,v1[2],v1[3],persp_res), "1", cex = 2)
text(trans3d(v2[1]-0.1,v2[2],v2[3],persp_res), "2", cex = 2)
text(trans3d(v3[1]-0.1,v3[2],v3[3],persp_res), "3", cex = 2)
text(trans3d(v4[1],v4[2],v4[3]+0.1,persp_res), "4", cex = 2)

###
# dev.copy(png,'myplot.png')
# dev.off()
dev.copy(pdf,"indepsurface.pdf")
dev.off()
