% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feasibility_check.R
\name{min_n_multistage}
\alias{min_n_multistage}
\title{Find the smallest total sample size such that a multi-stage trial meets desired constraints}
\usage{
min_n_multistage(args, cases, trial_method, objective_fun, min_n = 1,
  max_n = min_n * 1000, step_n = 10, showiter = FALSE)
}
\arguments{
\item{args}{a list containing a subset of the arguments for the functions \code{\link{getEffBounds}} and \code{\link{simTrial}} (or comparable functions, see \code{trial_method} argument). This should include a FWER constraint.}

\item{cases}{A list of power constraints, in the same format as those sent to \code{\link{optimizeTrial}}}

\item{trial_method}{either 'cov' or 'MB' for Maurer-Bretz (2013).}

\item{objective_fun}{see \code{\link{optimizeTrial}}}

\item{min_n}{The smallest sample size to consider}

\item{max_n}{The largest sample size to consider}

\item{step_n}{The step size to consider when carrying out the binary search. For example, if \code{step_n = 10}, \code{min_n=0}, and \code{max_n=100}, the function will find the smallest n_total satisfying the supplied constraints, and that is also a multiple of 10.}

\item{showiter}{passed to \code{\link{binsearch}}}
}
\value{
A list containing
\item{n}{The smallest feasible n_total}
\item{soln}{Output from \code{\link{binsearch}}}.
}
\description{
Takes as input a list of arguments (args) that define an adaptive trial (see \code{\link{buildTrial}} or \code{\link{optimizeTrial}}). This function adjusts the \code{n_total} argument in order to find the smallest maximum sample size that meets the power constraints specified in the \code{cases} argument.
}
\details{
This function requires that the objective function contain a 'base' element, and a 'power_diffs' element that is nonnegative when power constraints are met. For example, see \code{\link{min_E_SS_power_constraints}}.
}
\references{
Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
}
\seealso{
\code{\link{min_n_feasible}}, \code{\link{feasibility_check}}
}
