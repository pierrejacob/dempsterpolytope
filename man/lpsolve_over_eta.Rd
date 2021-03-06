% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpsolve_over_eta.R
\name{lpsolve_over_eta}
\alias{lpsolve_over_eta}
\title{Solve LP over polytope described by eta}
\usage{
lpsolve_over_eta(eta, objective_vector)
}
\arguments{
\item{eta}{A matrix of dimension K x K describing a convex polytope.}

\item{objective_vector}{a vector giving coefficient of linear function to minimize.}
}
\value{
Evaluation of the objective function at its minimum.
}
\description{
Minimize linear program of the form
t(objective_vector) * theta
over theta in the convex polytope of theta included in the simplex
and such that theta_l / theta_k <= eta\link{k,l}.
Employs the \code{lpSolveAPI} package.
}
\examples{
\dontrun{
eta <- rejection_sampler(c(1,0,3,2))$etas
lpsolve_over_eta(eta, c(1,0,0,0))
}
}
