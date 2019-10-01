#' ## function to run Coupled PIMH estimators
#' ## based on the output of the SMC sampler
#' ## nparticles: number of particles in the SMC sampler
#' ## X: sequence of observations
#' ## K: number of categories
#' ## resamplingtimes: resampling times for the SMC algorithm; can be obtained from a preliminary SMC run
#' ## k: equivalent of burnin for standard MCMC, in the sense that the unbiased estimators average from step k onwards
#' ## m: equivalent of the total number of iterations for standard MCMC, in the sense that unbiased estimators average from k to m
#' ## max_iterations: number of iterations after which to stop the while loop; default to Inf
#' ## h: test function of interest, default to NULL
#' ##
#' #'@export
#' coupled_pimh_chains <- function(nparticles, X, K, resamplingtimes, k = 0, m = 1, max_iterations = Inf, h = NULL){
#'   chain_state1 <- SMC_sampler(nparticles, X, K, resamplingtimes = resamplingtimes, h = h)
#'   chain_state2 <- SMC_sampler(nparticles, X, K, resamplingtimes = resamplingtimes, h = h)
#'   ## if test functions are provided, calculations need be performed
#'   mcmcestimator <- NULL
#'   biascorrection <- NULL
#'   dimh <- NULL
#'   if (!is.null(h)){
#'     mcmcestimator <- chain_state1$hestimator
#'     dimh <- length(mcmcestimator)
#'     if (k > 0){
#'       mcmcestimator <- rep(0, dimh)
#'     }
#'     biascorrection <- rep(0, dimh)
#'   }
#'   ##
#'   meetingtime <- Inf
#'   finished <- FALSE
#'   ### First step, propose chain_state2 in IMH for chain 1
#'   logu <- log(runif(1))
#'   imh_accept <- (logu < (sum(chain_state2$normcst) - sum(chain_state1$normcst)))
#'   if (imh_accept){
#'     meetingtime <- 1
#'     chain_state1 <- chain_state2
#'   } else {
#'     # do nothing
#'   }
#'   # at this point step t=1 is done, i.e. we have X_{t} and Y_{t-1}
#'   iter <- 1
#'   ## update estimator
#'   if (!is.null(h) && k <= 1){
#'     estimator1 <-
#'     if (k == 0){ ## then we need to compute h(X_1) - h(Y_0)
#'       biascorrection <- biascorrection + (min(1, (0 - k + 1)/(m - k + 1))) * (chain_state1$hestimator  - chain_state2$hestimator )
#'     }
#'     if (k <= iter && iter <= m){ ## we need to compute h(X_1) and update the MCMC estimator
#'       mcmcestimator <- mcmcestimator + chain_state1$hestimator
#'     }
#'   }
#'   # and we can now sample X_{t}, Y_{t-1} for t >= 2
#'   while (!finished && iter < max_iterations){
#'     # iter corresponds to t = 2,3,...
#'     iter <- iter + 1
#'     if (meetingtime <= iter){ # if already met, standard PIMH move
#'       smc_proposal <-   SMC_sampler(nparticles, X, K, resamplingtimes = resamplingtimes, h = h)
#'       logu <- log(runif(1))
#'       imh_accept <- (logu < (sum(smc_proposal$normcst) - sum(chain_state1$normcst)))
#'       if (imh_accept){
#'         chain_state1 <- smc_proposal
#'       }
#'       # update MCMC estimator, bias correction term is zero
#'       if ((!is.null(h)) && (k <= iter && iter <= m)){
#'         mcmcestimator <- mcmcestimator + chain_state1$hestimator
#'       }
#'     } else { # no meeting yet, so coupled PIMH move
#'       smc_proposal <-   SMC_sampler(nparticles, X, K, resamplingtimes = resamplingtimes, h = h)
#'       logu <- log(runif(1))
#'       imh_accept1 <- (logu < (sum(smc_proposal$normcst) - sum(chain_state1$normcst)))
#'       imh_accept2 <- (logu < (sum(smc_proposal$normcst) - sum(chain_state2$normcst)))
#'       if (imh_accept1){
#'         chain_state1 <- smc_proposal
#'       }
#'       if (imh_accept2){
#'         chain_state2 <- smc_proposal
#'       }
#'       if (imh_accept1 && imh_accept2){ # if both chains accept same proposal
#'         meetingtime <- iter
#'       }
#'       if (!is.null(h)){ # update MCMC estimator and bias correction
#'         if ((k <= iter) && (iter <= m)){ # update sum_{s=k}^t h(X_s) with last element h(X_t)
#'           mcmcestimator <- mcmcestimator + chain_state1$hestimator
#'         }
#'         if ((k + 1) <= iter){
#'           # update sum_{s=k+1}^t [h(X_{s}) - h(Y_{s-1})]
#'           biascorrection <- biascorrection + (min(1, (iter - k)/(m - k + 1))) * (chain_state1$hestimator - chain_state2$hestimator)
#'         }
#'       }
#'     }
#'     # stop after max(m, tau) steps
#'     if (iter >= max(meetingtime, m)){
#'       finished <- TRUE
#'     }
#'   }
#'   mcmcestimator <- mcmcestimator / (m - k + 1)
#'   uestimator <- mcmcestimator + biascorrection
#'   return(list(mcmcestimator = mcmcestimator, biascorrection = biascorrection, uestimator = uestimator,
#'               meetingtime = meetingtime, iteration = iter, finished = finished))
#' }
