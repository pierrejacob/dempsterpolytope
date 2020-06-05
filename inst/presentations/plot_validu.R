## This script generates a plot illustrating how to construct a valid 'u' in R_x
rm(list = ls())
library(dempsterpolytope)
library(latex2exp)
graphsettings <- set_custom_theme()
set.seed(1)

attach(graphsettings)

g <- create_plot_triangle(graphsettings)
pt_bar <- c(0.3, 0.25, 0.45)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)
v1 <- v_cartesian[[1]]
v2 <- v_cartesian[[2]]
v3 <- v_cartesian[[3]]

g <- create_plot_triangle(graphsettings)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 1, fill = graphsettings$cols[1], alpha = 0.7)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 2, fill = graphsettings$cols[2], alpha = 0.7)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 3, fill = graphsettings$cols[3], alpha = 0.7)
meanpi1 <- c(pt_xy[1]/3 + v2[1]/3 + v3[1]/3, pt_xy[2]/3 + v2[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(Delta[1](theta)), parse = FALSE, size = 5) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(Delta[2](theta)), parse = FALSE, size = 5) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(Delta[3](theta)), parse = FALSE, size = 5) 
g <- g + annotate(geom = 'text', size = 5, x = .56, y = .28, label = TeX("$\\theta$", output = "character"), parse = TRUE)

K <- 3

df <- data.frame()
counts <- c(2,3,1)
for (d in 1:K){
  pts_barcoord <- dempsterpolytope:::runif_piktheta_cpp(counts[d], d, pt_bar)$pts
  pts_cart <- t(apply(pts_barcoord, 1, function(v) barycentric2cartesian(v, v_cartesian)))
  df <- rbind(df, data.frame(x = pts_cart[,1], y = pts_cart[,2], category = d))
}
df
gpoints <- g + geom_point(data=df, aes(x = x, y = y, fill = factor(category), col = factor(category)), size = 3, shape = 21)
gpoints <- gpoints + scale_fill_manual("", values = cols) + scale_color_manual("", values = contcols) + theme(legend.position = "none")
gpoints


ggsave(filename = "validu.pdf", plot = gpoints, width = 5, height = 5)

