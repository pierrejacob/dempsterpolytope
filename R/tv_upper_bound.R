## from vector of meeting times, lag L and time t, obtain upper bound on TV
#'@export
tv_upper_bound <- function(coupling_times, L, t){
  return(mean(pmax(0,ceiling((coupling_times-L-t)/L))))
}
