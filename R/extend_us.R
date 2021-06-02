## takes "Us" as in the output of 'gibbs_sampler'
## and adds components so as to add empty categories
## Us should be a list with d entries
## each being an array of size 'niterations' x 'N_k' x 'd'
## where N_k corresponds to the counts for this category, d the number of non-empty categories
## whichbefore and whichnew should be vector of unique indices covering 1,...,K, where K 
## is the number of categories including the empty ones
#'@export
extend_us <- function(Us, whichbefore, whichnew){
  # retrieve number of iterations
  niterations <- dim(Us[[1]])[1]
  # retrieve counts from existing Us
  countsbefore <- unlist(sapply(Us, function(l){
    count_k <- dim(l)[2]
    if (is.null(count_k)) return(0)
    else return(count_k)
  }
  ))
  ## *important note:
  ## *at this point the function assumes these counts are all non zero
  ## *could be changed later...
  stopifnot(length(countsbefore) == length(whichbefore))
  ## total number of categories
  K <- max(c(whichbefore, whichnew))
  counts <- rep(0, K)
  counts[whichbefore] <- countsbefore
  ## c(whichbefore, whichnew) should be a partition of 1:K
  stopifnot(length(unique(c(whichbefore, whichnew))) == length(c(whichbefore, whichnew)))
  Kbefore <- length(whichbefore)
  ## now we need to add back the empty categories
  ## function takes a matrix of Us and extends them 
  extend_u <- function(us){
    if (is.null(dim(us))) us <- matrix(us, nrow = 1)
    s <- rgamma(dim(us)[1], shape = Kbefore, rate = 1)
    expos <- matrix(rexp(dim(us)[1] * length(whichnew), rate = 1), nrow = dim(us)[1])
    extendedu <- matrix(0, nrow = dim(us)[1], ncol = K)
    extendedu[,whichbefore] <- us * s
    extendedu[,whichnew] <- expos
    extendedu <- extendedu / (s + rowSums(expos))
    return(extendedu)
  }
  ## next extend all the U's
  ## and compute associated etas
  extendedUs <- list()
  extendedetas <- array(Inf, dim = c(niterations, K, K))
  inonempty <- 0
  for (k in 1:K){
    if (counts[k] == 0){
      extendedUs[[k]] <- NA
    } else {
      ## where was the Us corresponding to this k before extension
      indexbefore <- order(c(whichbefore, whichnew))[k]
      # inonempty <- inonempty + 1
      extendedUs[[k]] <- array(0, dim = c(niterations, counts[k], K))
      for (iteration in 1:niterations){
        ## extend Us
        extendedUs[[k]][iteration,,] <- extend_u(Us[[indexbefore]][iteration,,])
        ## recompute etas from Us
        extendedetas[iteration,k,] <- sapply(1:K, function(j) {
          return(min(extendedUs[[k]][iteration,,j]/extendedUs[[k]][iteration,,k]))
        })
      }
    }
  }
  return(list(Us = extendedUs, etas = extendedetas))
}
