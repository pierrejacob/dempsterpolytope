## This script generates the two plots of Figure 1 of the article.
library(dempsterpolytope)
library(latex2exp)

set_custom_theme()
set.seed(1)
rm(list = ls())

## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "yellow", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
triangle.df <- data.frame(x = c(v1[1], v2[1], v3[1]), y = c(v1[2], v2[2], v3[2]))
##
K <- 3
matrixT <- matrix(0, nrow = K-1, ncol = K-1)
for (k in 1:(K-1)){
  kth_components <- sapply(v_cartesian, function(x) x[k])
  matrixT[k,] <- kth_components[-K] - kth_components[K]
}
#

pt_bar <- c(0.3, 0.25, 0.45)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)
# Now with ggplot2
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
##

g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v2[1], v3[1]), y = c(pt_xy[2], v2[2], v3[2])),
                      aes(fill = "1"), colour = "black", alpha = 0.5)
g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v3[1]), y = c(pt_xy[2], v1[2], v3[2])),
                      aes(fill = "2"), colour = "black", alpha = 0.5)
g <- g + geom_polygon(data = data.frame(x = c(pt_xy[1], v1[1], v2[1]), y = c(pt_xy[2], v1[2], v2[2])),
                      aes(fill = "3"), colour = "black", alpha = 0.5)
g <- g + scale_fill_manual(name = "partition: ", values = cols) + theme(legend.position = "none")
meanpi1 <- c(pt_xy[1]/3 + v2[1]/3 + v3[1]/3, pt_xy[2]/3 + v2[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(Delta[1](theta)), parse = FALSE, size = 5) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(Delta[2](theta)), parse = FALSE, size = 5) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(Delta[3](theta)), parse = FALSE, size = 5) 
g <- g + annotate(geom = 'text', size = 5, x = .56, y = .28, label = TeX("$\\theta$", output = "character"), parse = TRUE)
g
# ggsave(filename = "subsimplices.pdf", plot = g, width = 5, height = 5)







