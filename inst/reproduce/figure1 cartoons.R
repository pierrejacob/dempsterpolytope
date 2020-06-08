## This script generates the two plots of Figure 1 of the article.
rm(list = ls())
library(dempsterpolytope)
library(latex2exp)
graphsettings <- set_custom_theme()
set.seed(1)

# number of categories
K <- 3
categories <- 1:K
# data 
counts <- c(2,3,1)

## triangle with equal sides
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
pt_bar <- c(0.3, 0.25, 0.45)
pt_xy  <- barycentric2cartesian(pt_bar, v_cartesian)
# Now with ggplot2
g <- create_plot_triangle(graphsettings)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 1, fill = graphsettings$cols[1], alpha = 0.8)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 2, fill = graphsettings$cols[2], alpha = 0.8)
g <- add_plot_subsimplex(graphsettings, g, pt_bar, 3, fill = graphsettings$cols[3], alpha = 0.8)
meanpi1 <- c(pt_xy[1]/3 + v2[1]/3 + v3[1]/3, pt_xy[2]/3 + v2[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi1[1], y = meanpi1[2], label = expression(Delta[1](theta)), parse = FALSE, size = 7) 
meanpi2 <- c(pt_xy[1]/3 + v1[1]/3 + v3[1]/3, pt_xy[2]/3 + v1[2]/3 + v3[2]/3)
g <- g + annotate(geom = "text", x = meanpi2[1], y = meanpi2[2], label = expression(Delta[2](theta)), parse = FALSE, size = 7) 
meanpi3 <- c(pt_xy[1]/3 + v2[1]/3 + v1[1]/3, pt_xy[2]/3 + v2[2]/3 + v1[2]/3)
g <- g + annotate(geom = "text", x = meanpi3[1], y = meanpi3[2], label = expression(Delta[3](theta)), parse = FALSE, size = 7) 
g <- g + annotate(geom = 'text', size = 7, x = .56, y = .28, label = TeX("$\\theta$", output = "character"), parse = TRUE)
g
ggsave(filename = "subsimplices.pdf", plot = g, width = 5, height = 5)




## triangle with equal sides
attach(graphsettings)
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
# Now with ggplot2
g <- create_plot_triangle(graphsettings)

pt_bar1 <- c(0.3, 0.25, 0.45)
g <- add_plot_subsimplex(graphsettings, g, pt_bar1, 1, direction = "reverse", fill = graphsettings$cols[1], 
                         colour = graphsettings$contcols[1], alpha = 0.8)
g <- add_plot_points(graphsettings, g, pt_bar1)

pt_bar2 <- c(0.4, 0.5, 0.1)
g <- add_plot_subsimplex(graphsettings, g, pt_bar2, 2, direction = "reverse", fill = graphsettings$cols[2], 
                         colour = graphsettings$contcols[2], alpha = 0.8)
g <- add_plot_points(graphsettings, g, pt_bar2)

pt_bar3 <- c(0.1, 0.45, 0.45)
g <- add_plot_subsimplex(graphsettings, g, pt_bar3, 3, direction = "reverse", fill = graphsettings$cols[3], 
                         colour = graphsettings$contcols[3], alpha = 0.8)
g <- add_plot_points(graphsettings, g, pt_bar3)

g <- g + annotate(geom = 'text', size = 7, x = .6, y = .2, label = TeX("$u_1$", output = "character"), parse = TRUE)
g <- g + annotate(geom = 'text', size = 7, x = .45, y = .1, label = TeX("$u_2$", output = "character"), parse = TRUE)
g <- g + annotate(geom = 'text', size = 7, x = .3, y = .4, label = TeX("$u_3$", output = "character"), parse = TRUE)
g

ggsave(filename = "cups.pdf", plot = g, width = 5, height = 5)


