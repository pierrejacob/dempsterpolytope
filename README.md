
## dempsterpolytope

Implements a Gibbs sampler for Dempster’s inference approach for
Categorical distributions; this package is a companion to an article by
Pierre E. Jacob, Ruobin Gong, Paul T. Edlefsen, Arthur P. Dempster. The
link to an online version of the article will be added here at some
point.

This is not a general-purpose statistical software. This is just a
collection of scripts intended to reproduce figures and tables of a
paper. Use at your own risk\!

The folder inst/reproduce/ contains the scripts to reproduce the
figures. The folder inst/tests/ contains internal checks, and vignettes/
contains tutorials on how to use the package’s main functions, which at
the moment are just R scripts and not proper “R vignettes”.

### Installation

The package can be installed from R via:

``` r
# install.packages("devtools")
devtools::install_github("pierrejacob/dempsterpolytope")
```

It depends on the packages Rcpp, RcppEigen, igraph, rcdd, dplyr,
lpSolveAPI which can be installed
via:

``` r
install.packages(c("Rcpp", "RcppEigen", "igraph", "rcdd", "dplyr", "lpSolveAPI"))
```

Additionally you might want to install other packages, to help with
parallel computation:

``` r
install.packages(c("doParallel", "doRNG"))
```

and to help with manipulating results and plotting:

``` r
install.packages(c("tidyr", "ggplot2"))
```

although these packages are not strictly required.

### Usage

The “Dempster-Shafer” approach to inference leads, in the case of
Categorical distributions, to random convex polytopes. The Gibbs sampler
implemented in this package generates such polytopes. The following code
shows a random polytope obtained after some iterations of the Gibbs
sampler, for the data set (12, 4, 7), meaning 12 observations in
category 1, 14 in category 2, 7 in category 3.

``` r
library(dempsterpolytope)
set.seed(1)

# count data
counts <- c(10, 4, 7)
# number of MCMC iterations
niterations <- 100
# run Gibbs sampler
gibbs_results <- gibbs_sampler(niterations, counts)

# obtain a K x K matrix representing a convex polytope in the simplex
# by taking the terminal iteration of the Gibbs chain
eta <- gibbs_results$etas_chain[niterations,,]
# convert polytope to H-representation and V-representation
eta_converted <- etas2cvxpolytope(eta)
# H-representation, or "half-plane" representation
# means the set is represented as points x such that 'constr' times x <= 'rhs' 
eta_converted$constr
#> $constr
#>            [,1]      [,2]
#>  [1,]  1.000000  1.000000
#>  [2,] -1.000000  0.000000
#>  [3,]  0.000000 -1.000000
#>  [4,] -0.344655  1.000000
#>  [5,] -1.608173 -1.000000
#>  [6,]  1.000000 -3.351633
#>  [7,] -1.000000 -6.775739
#>  [8,]  3.394207  2.394207
#>  [9,]  1.271098  2.271098
#> 
#> $rhs
#> [1]  1.000000  0.000000  0.000000  0.000000 -1.000000  0.000000 -1.000000
#> [8]  2.394207  1.271098
#> 
#> $dir
#> [1] "<=" "<=" "<=" "<=" "<=" "<=" "<=" "<=" "<="
# V-representation, or "vertex" representation
# means the set is represented by its vertices
eta_converted$vertices_barcoord
#>           [,1]      [,2]      [,3]
#> [1,] 0.5120778 0.1764902 0.3114321
#> [2,] 0.5245116 0.1564944 0.3189940
#> [3,] 0.5827381 0.1738669 0.2433950
#> [4,] 0.5674307 0.1955678 0.2370015

# next we can view the polytope in the K-simplex as a triangle, since here K = 3
v_cartesian <- list(c(1/2, sin(pi/3)), c(0,0), c(1,0))
cols <- c("red", "green", "blue")
v1 <- v_cartesian[[1]]; v2 <- v_cartesian[[2]]; v3 <- v_cartesian[[3]]
## this loads ggplot2
set_my_theme()
g <- ggplot_triangle(v_cartesian, etas = eta, addpolytope = T, cols = cols)
g
```

![](README-usage-1.png)<!-- -->

The generated polytope is the black polygon in the middle of the
triangle, which represents the entire simplex.
