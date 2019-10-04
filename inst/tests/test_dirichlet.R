## This script tests calculations to do with Dirichlet distributions

library(dempsterpolytope)
set_my_theme()
## triangle with equal sides
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]

## going from K = 2 to K = 3

npts <- 10000
## 
u1 <- runif(npts, 0, 1)
u2 <- 1 - u1
s <- rgamma(npts, 2, 1)
w <- rexp(npts, 1)
samples2 <- cbind(u1, u2, s, w)
# plot(x = samples2[,1], y = samples2[,2])

## transformation
transf <- function(u1, u2, s, w){
  return(c(s * u1 / (s + w), s * u2 / (s + w), w / (s + w), s + w))
}
invtransf <- function(z1, z2, z3, z4){
  return(c(z1 / (1 - z3), z2 / (1 - z3), (1-z3) * z4, z3 * z4))
}
samples2transf <- t(apply(samples2, 1, function(v) transf(v[1], v[2], v[3], v[4])))
samples2untransf <- t(apply(samples2transf, 1, function(v) invtransf(v[1], v[2], v[3], v[4])))
head(samples2 - samples2untransf)
## ok so inv transformation is valid
## Now are the first three components of samples2transf distributed according to a Dirichlet ?
# 1st and 2nd components
plot(x = samples2transf[,1], y = samples2transf[,2], col = rgb(0,0,0,0.1))
hist(samples2transf[,1], nclass = 100, prob = TRUE)
curve(dbeta(x, 1, 2), add = TRUE)
hist(samples2transf[,2], nclass = 100, prob = TRUE)
curve(dbeta(x, 1, 2), add = TRUE)
# 4th component
hist(samples2transf[,4], nclass = 100, prob = TRUE)
curve(dgamma(x, 3, 1), add = TRUE)
## it looks like it's correct...


##
K <- 3
rdiri <- function(K){
  x <- rexp(K)
  return(x / sum(x))
}
npts <- 10000
pts_barcoord <- matrix(nrow = npts, ncol = K)
for (n in 1:npts){
  pts_barcoord[n,] <- rdiri(K)
}

pts_cartesian <- t(apply(pts_barcoord, 1, function(x) barycentric2cartesian(x, v_cartesian)))
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
g <- g + geom_point(data=data.frame(x = pts_cartesian[,1], y = pts_cartesian[,2]), alpha = 0.25)
g

hist(pts_barcoord[,1], nclass = 100, prob = TRUE)
curve(dbeta(x, 1, 2), add = TRUE)

### go to Dirichlet in smaller dimension

pts_barcoord_smaller <- pts_barcoord[,1:2]
pts_barcoord_smaller <- t(apply(pts_barcoord_smaller, 1, function(v) v / sum(v)))
hist(pts_barcoord_smaller[,1], nclass = 100, prob = TRUE)
hist(pts_barcoord_smaller[,2], nclass = 100, prob = TRUE)

## go back to Dirichlet in larger dimension?

## sum of K exponenient
# hist(sapply(1:1e4, function(x) sum(rexp(K))), prob = TRUE, nclass = 100)
# curve(dgamma(x, K, 1), add = TRUE)
sums_ <- rgamma(nrow(pts_barcoord_smaller), 2, 1)
y1 <- pts_barcoord_smaller[,1] * sums_
y2 <- pts_barcoord_smaller[,2] * sums_
y3 <- rexp(nrow(pts_barcoord_smaller), 1)
pts_barcoord_larger2 <- cbind(y1, y2, y3)
pts_barcoord_larger2 <- t(apply(pts_barcoord_larger2, 1, function(v) v / sum(v)))
pts_cartesian2 <- t(apply(pts_barcoord_larger2, 1, function(x) barycentric2cartesian(x, v_cartesian)))

##
g <- ggplot(triangle.df, aes(x = x, y = y)) + geom_polygon(fill = "white", colour = "black")
g <- g + scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + xlab("") + ylab("")
g <- g + geom_text(data = data.frame(x = v2[1], y = v2[2]-0.1, label = "2"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v3[1], y = v3[2]-0.1, label = "3"), aes(label = label), size = 5)
g <- g + geom_text(data = data.frame(x = v1[1], y = v1[2]+0.1, label = "1"), aes(label = label), size = 5)
g <- g + geom_point(data=data.frame(x = pts_cartesian2[,1], y = pts_cartesian2[,2]), alpha = 0.25)
g

hist(pts_barcoord_larger2[,1], nclass = 100, prob = TRUE)
curve(dbeta(x, 1, 2), add = TRUE)

hist(pts_barcoord_larger2[,2], nclass = 100, prob = TRUE)
curve(dbeta(x, 1, 2), add = TRUE)


