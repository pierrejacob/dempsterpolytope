## minimize t(objective_vector) * theta over theta in simplex such that theta_l / theta_k <= eta[k,l]
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
  ##
  ## could set the lower bounds to be -Inf instead of the default 0, but no difference here
  ## lpSolveAPI::set.bounds(lpobject, lower = rep(-Inf, ncol(lc$constr)))
  ##
  ## find solution
  solve(lpobject)
  ## save objective function evaluated at the solution
  result <- lpSolveAPI::get.objective(lpobject)
  ## remove LP object (not sure why)
  rm(lpobject)
  return(result)
}
