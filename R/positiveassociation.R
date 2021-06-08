#'@export
positiveassociation <- function(eta){
  ## check relation between feasible set, i.e. polyope of theta s.t. theta_j / theta_d <= etas[d,j] 
  ## with positive association assumption
  ## -log(theta_1) - log(theta_4) + log(theta_2) + log(theta_3)  <= 0
  ## this function returns two booleans
  ## 1) whether feasible set intersects with the 'positive association polytope' 
  ## 2) whether feasible set is contained in 'positive association polytope'
  if (dim(eta)[1] != 4){
    stop("The matrix etas must be 4x4 for the function 'positiveassociation'")
  }
  baryvertices <- etas_vertices(eta)
  ## round to zero if close to zero
  baryvertices[abs(baryvertices)<1e-10] <- 0
  ## compute 
  positivassoc_transformation <- apply(baryvertices, 1, function(v) -log(v[1]) -log(v[4])+log(v[2])+log(v[3]))
  positivassoc_transformation[is.na(positivassoc_transformation)] <- +Inf
  return(c(any(positivassoc_transformation<0), all(positivassoc_transformation<0)))
}

