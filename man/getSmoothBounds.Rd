% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildTrial.R
\name{getSmoothBounds}
\alias{getSmoothBounds}
\title{Get smooth parametric boundaries for efficacy and futility}
\usage{
getSmoothBounds(n_per_stage, delta_vec, coefficients_vec, intercepts_vec,
  takeDiffs)
}
\arguments{
\item{n_per_stage}{vector of sample sizes, or a vector proportional to these sample sizes}

\item{delta_vec}{a vector with elements greater than or equal to 0, one for each boundary to be computed (e.g. one for \eqn{H01}, one for \eqn{H02}, and one for \eqn{H0C}).}

\item{coefficients_vec}{a vector of constants, one for each boundary to be computed (e.g. one for \eqn{H01}, one for \eqn{H02}, and one for \eqn{H0C})}

\item{intercepts_vec}{a vector of constants to be added to the set to zero for efficacy boundaries}

\item{takeDiffs}{set to TRUE if calculating efficacy boundaries}
}
\value{
A matrix with each row corresponding to one boundary (one hypothesis).
}
\description{
A general parametric form for allocating alpha or creating futility boundaries. For details, see Fisher and Rosenblum (2016).
}
\examples{
getSmoothBounds(
	  n_per_stage=1:5,
	  intercepts_vec=c(0,0,1,1,0,1),
	  delta_vec=c(1,1,1,1,1,1),
	  coefficients_vec=c(1,1,1,1,2,2),
  takeDiffs=FALSE
	  )
}
\references{
Aaron Fisher and Michael Rosenblum (2016). Stochastic Optimization of Adaptive Enrichment Designs for Two Subpopulations. http://biostats.bepress.com/jhubiostat/paper279/
}
