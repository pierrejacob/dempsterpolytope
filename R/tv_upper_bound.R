## implements the TV upper bounds of https://arxiv.org/abs/1905.09971
## based on meeting times
## from vector of meeting times, lag L and time t, obtain upper bound on TV
#'@export
tv_upper_bound <- function(meetingtimes, L, t){
  return(mean(pmax(0,ceiling((meetingtimes-L-t)/L))))
}
