# for a given phi_0, and assertion {phi < phi_0}
# upper probability is proportion of time random interval [phi_lower, phi_upper] intersects with [0, phi_0]
# lower probability is proportion of time random interval [phi_lower, phi_upper] is contained in [0, phi_0]
#'@export
linkage_cdf_lowerupper <- function(phi_0, lu_chain){
  # intersect if phi_lower < phi_0
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < phi_0))
  # contains if phi_upper < phi_0
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < phi_0))
  return(c(cdf_lower, cdf_upper))
}

