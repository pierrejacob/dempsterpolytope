#'@rdname lpsolve_over_eta
#'@title Solve LP over polytope described by eta
#'@description Minimize linear program of the form
#'t(objective_vector) * theta 
#'over theta in the convex polytope of theta included in the simplex
#'and such that theta_l / theta_k <= eta[k,l]. 
#'Employs the \code{lpSolveAPI} package.
#'@param eta A matrix of dimension K x K describing a convex polytope.
#'@param objective_vector a vector giving coefficient of linear function to minimize.
#'@return Evaluation of the objective function at its minimum.
#'@examples
#'\dontrun{
#'eta <- rejection_sampler(c(1,0,3,2))$etas
#'lpsolve_over_eta(eta, c(1,0,0,0))
#'}
#'@export
lpsolve_over_eta <- function(eta, objective_vector){
  ## get constraints from simplex
  lc1 <- simplex_linearconstraints(dim(eta)[1])
  ## get constraints from eta, i.e. theta_l / theta_k <= eta[k,l]
  lc2 <- eta_linearconstraints(eta)
  ## combine constraints
  lc <- concatenate_linearconstraints(lc1, lc2)
  ## create linear program 'LP' object
  lpobject <- lpSolveAPI::make.lp(nrow = dim(lc$constr)[1], 
                                  ncol = dim(lc$constr)[2])
  ## set constraints
  lpSolveAPI::set.constr.type(lpobject, lc$dir)
  lpSolveAPI::set.rhs(lpobject, lc$rhs)
  for (icol in 1:ncol(lc$constr)){
    lpSolveAPI::set.column(lpobject, icol, lc$constr[,icol])
  }
  ## set objective function 't(objective_vector) theta'
  lpSolveAPI::set.objfn(lpobject, objective_vector)
  ## find solution
  solve(lpobject)
  ## save objective function evaluated at the solution
  result <- lpSolveAPI::get.objective(lpobject)
  ## remove LP object (not sure why)
  rm(lpobject)
  return(result)
}

#'@rdname lpsolve_over_eta_log
#'@title Solve LP over polytope described by log(eta)
#'@description Minimize linear program of the form
#'t(objective_vector) * log(theta) 
#'over theta in the convex polytope of log(theta) 
#' such that log(theta_l) - log(theta_k) <= log(eta[k,l])
#' and such that sum_k log (theta_k) = 0. 
#'Employs the \code{lpSolveAPI} package.
#'@param eta A matrix of dimension K x K describing a convex polytope.
#'@param objective_vector a vector giving coefficient of linear function (of log theta) to minimize.
#'@return Evaluation of the objective function at its minimum.
#'@examples
#'\dontrun{
#'eta <- rejection_sampler(c(1,0,3,2))$etas
#'lpsolve_over_eta_log(eta, c(1,0,0,0))
#'}
#'@export
lpsolve_over_eta_log <- function(eta, objective_vector){
  constraints <- eta_log_linearconstraints(eta)
  ## create linear program 'LP' object
  lpobject <- lpSolveAPI::make.lp(nrow = dim(constraints$constr)[1], 
                                  ncol = dim(constraints$constr)[2])
  ## set constraints
  lpSolveAPI::set.constr.type(lpobject, constraints$dir)
  lpSolveAPI::set.rhs(lpobject, constraints$rhs)
  for (icol in 1:ncol(constraints$constr)){
    lpSolveAPI::set.column(lpobject, icol, constraints$constr[,icol])
  }
  ## set objective function 't(objective_vector) log(theta)'
  lpSolveAPI::set.objfn(lpobject, objective_vector)
  lpSolveAPI::set.bounds(lpobject, lower = rep(-Inf, ncol(constraints$constr)))
  ## find solution
  solve(lpobject)
  ## save objective function evaluated at the solution
  result <- lpSolveAPI::get.objective(lpobject)
  ## remove LP object (not sure why)
  rm(lpobject)
  return(result)
}

