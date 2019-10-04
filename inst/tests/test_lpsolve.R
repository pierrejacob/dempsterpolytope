## This scripts tests that the solution of a linear program
## matches the solution of a shortest path problem 
library(dempsterpolytope)
set.seed(1)
rm(list = ls())
library(lpSolveAPI)

## number of categories
K <- 3
## number of observations
n <- 5
freqX <- 1:K

## run Gibbs sampler to get random polytopes
niterations_gibbs <- 1e2
samples_gibbs <- gibbs_sampler(niterations_gibbs, freqX)
etas <- samples_gibbs$etas_chain[44,,]

# now suppose we want to optimize theta_k
# under the constraint that theta is in the simplex
# and that theta_i / theta_j < eta_{j,i} for all i, and for j not equal to k

kupdate <- 2

# the program looks likes
# minimize      - x_k
# subject to      x_1 + x_2 + x_3   = 1
#                 x_1              >= 0
#                       x_2        >= 0
#                             x_3  >= 0
#                 

# number of constraints = K+1
# number of decision variables = 3
my.lp <- make.lp(K+1, K)
vec_ <- c(1, rep(0, K))
for (i in 1:K){
  vec <- vec_
  vec[i+1] <- 1
  set.column(my.lp, i, vec)
}

vec_ <- rep(0, K)
vec_[kupdate] <- -1
set.objfn(my.lp, vec_)
set.rhs(my.lp, c(1, rep(0, K)))
set.constr.type(my.lp, c("=", rep(">=", K)))
#
## add i.e. theta_i - eta_{j,i} theta_j < 0 
for (j in setdiff(1:K, kupdate)){
  for (i in setdiff(1:K, j)){
    cst <- rep(0, K)
    cst[i] <- 1
    cst[j] <- -etas[j,i]
    add.constraint(my.lp, cst, "<=", 0)
  }
}
my.lp
solve(my.lp)
# get.objective(my.lp)
theta_star <- get.variables(my.lp)
# get.constraints(my.lp)

## compare to graph-based approach
g <- graph_from_adjacency_matrix(log(etas), mode = "directed", weighted = TRUE, diag = FALSE)
theta_star_graph <- rep(0, K)
notk <- setdiff(1:K, kupdate)
minimum_values <- rep(1, K)
minimum_values[notk] <- distances(g, v = notk, to = kupdate, mode = "out")[,1]
theta_star_graph <- exp(-minimum_values)
theta_star_graph[kupdate] <- 1
theta_star_graph <- theta_star_graph / sum(theta_star_graph)
print(theta_star_graph)
# versus
print(theta_star)

## define the same LP but using "set.column" only, as the "add.constraint" should be
## avoided according to the help page of "add.constraint".

## number of constraints:
# 1 for sum of component equals 1
# K for each component >= 0
# then for the update given category k, K-1 other categories give K-1 constraints each
# so in total 1 + K + (K-1)*(K-1)
Km1squared <- (K-1)*(K-1)
nconstraints <- K + 1 + Km1squared
mat_cst <- matrix(0, nrow = nconstraints, ncol = K)
mat_cst[1,] <- 1
for (i in 1:K) mat_cst[1+i,i] <- 1
mat_cst
dir <- c("=", rep(">=", K), rep("<=", Km1squared))
rhs <- c(1, rep(0, K), rep(0, Km1squared))
# then for a given etas, if we are updating component kupdate
icst <- 1
for (j in setdiff(1:K, kupdate)){
  for (i in setdiff(1:K, j)){
    # theta_i - eta_{j,i} theta_j < 0 
    row_ <- (K+1)+icst
    mat_cst[row_,i] <- 1
    mat_cst[row_,j] <- -etas[j,i]
    icst <- icst + 1
  }
}

#
directlp <- make.lp(nrow = nconstraints, ncol = K)
for (i in 1:K){
  set.column(directlp, i, mat_cst[,i])
}
vec_ <- rep(0, K)
vec_[kupdate] <- -1
set.objfn(directlp, vec_)
set.rhs(directlp, rhs)
set.constr.type(directlp, dir)
solve(directlp)
theta_star2 <- get.variables(directlp)

theta_star
theta_star2
theta_star_graph



