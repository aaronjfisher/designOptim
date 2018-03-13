# Functions for checking feasibility of a design,
# and searching for the closest design that is feasible.





###### Helper functions
######

# Interleave two vectors, v1 and v2.
# [INTERNAL FUNCTION]
# This code is based on code from Bogdan Romocea, at https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
# @param v1 one vector to interleave
# @param v2 a second vector to interleave
interleave <- function(v1,v2){
  if(length(v1)==0) return(v2)
  if(length(v2)==0) return(v1)

  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

# Reorder a vector so that mid points come first. This is similar in spirit to binary search.
# This recursive function is called exactly as many times as the length of \code{x}.
# [INTERNAL FUNCTION]
# @param x is a vector of integers returned by order(x)
orderByMids<-function(x){

  ox <- order(x)

  move <- ceiling(length(x)/2)
  xm<- x[move]
  x1<- x[ox<move]
  x2<- x[ox>move]

  if(length(x1) ==0 & length(x2)==0) return(xm)

  if(length(x1)>0) x1<-orderByMids(x1)
  if(length(x2)>0) x2<-orderByMids(x2)
 
  return(c(xm,interleave(x1,x2)))

}

######
######




#' Check feasibility of power constraints
#'
#' Functions to check feasibility of a design, adjust the design in order to make it feasible, or adjust the constraints so that they can be achieved. Here, "feasibility" means that there exists a 1-stage design at the supplied sample size that meets the supplied power and Type I error constraints.
#'
#' \code{feasibility_check} checks whether it is possible to meet a set of power constraints while maintaining a given familywise Type I error rate (FWER) and total sample size.
#'
#' \code{min_n_feasible} implements a binary search to find the smallest sample size that meets the power constraints and FWER specified. In each iteration of the search, \code{min_n_feasible} calls \code{feasibility_check}.
#'
#' Given a sample size and set of power cases, \code{max_power_feasible} finds the number \eqn{m} between 0 and 1 such that if all of the minimum power thresholds are multiplied by \eqn{m}, then these power constraints will be satisfied at the supplied value of \code{n_total}. In other words, it finds the factor by which the power constraints must be relaxed.
#'
#' @param FWER the required familywise Type I error rate for the trial
#' @param p1 population proportion in subpopulation 1
#' @param trial_method the type of trial to run. 'cov' for covariance based, 'MB' for Maurer Bretz, and 'covMB' for a combination approach.
#' @param r1 probability of being randomized to treatment in subpopulation 1
#' @param r2 probability of being randomized to treatment in subpopulation 2
#' @param n_total the total sample size for the 1-stage trial
#' @param cases a list of power constraints, of the same format as those passed to \code{\link{optimizeTrial}}
#' @param npoints_sqrt \code{feasibility_check} determines feasibility by searching over a grid of points for the alpha to be allocated between the two trials. Determines the number of grid points to search over. The search will be conducted over a triangle of points, with \code{(npoints_sqrt^2)/2} points.
#' @param trial_args either an empty list (for trial_method=='cov'), or a list containing graph edges for alpha reallocation (for trial_method=='MB').
#' @return \code{feasibility_check} returns a named vector, with the first element being 1 if the trial is feasible, and zero otherwise. The remaining elements of this vector tell the alpha allocated to each hypothesis in one particular setup that results in a feasible trial. \cr \code{min_n_feasible} returns the smallest sample size that meets the power constraints and FWER specified, as well as the output from \code{feasibility_check} at that sample size. \cr \code{max_power_feasible} returns the multiplier \eqn{m} and a list of modified cases in which each power threshold has been multiplied by \eqn{m}.
#' @export
#'
feasibility_check <- function(FWER, p1, trial_method, r1=0.5, r2=1-r1, n_total,cases,npoints_sqrt=10,trial_args=list()){

seq1 <- seq(10^-10,FWER,length=npoints_sqrt)
seq1 <- seq1[orderByMids(order(seq1))]
  #reordering the grid search in a smart way can let us stop the search sooner if we find a feasible design.
if(npoints_sqrt==1) seq1 <- FWER / 3 #equal allocation to each hypothesis

for(alpha_1 in seq1){
  seq2<- seq1[seq1<=FWER-alpha_1]

  for(alpha_2 in seq2){
  
    # Reset `all_constraints_satisfied` to TRUE
    all_constraints_satisfied <- TRUE
 
    for(case in cases){
      # Set means and variances according to user input for power constraint:
        #!! In future we can make this code cleaner by combining many of these into trial_args
        # (specifically p1, r1, r2, FWER), but for now we're keeping it just to avoid 
        # compatibility errors.

      trial_args$iter <- 10000 #!! Hard-coded
      trial_args$FWER <- FWER
      trial_args$num_stages<-1
      trial_args$n_total<-n_total
      trial_args$n_per_stage<-c(n_total)
      trial_args$p1<-p1
      trial_args$r1<-r1
      trial_args$r2<-r2
      trial_args$H01_eff_allocated<-alpha_1
      trial_args$H02_eff_allocated<-alpha_2
      trial_args$H0C_eff_allocated<-(FWER-(alpha_1+alpha_2))
      trial_args$H01_futility_boundaries<- -Inf # these are essentially ignored anyway
      trial_args$H02_futility_boundaries<- -Inf
      trial_args$H0C_futility_boundaries<- -Inf
      trial_args$time_limit<-200
      trial_args$enrollment_rate_combined<-100 
      trial_args$delay<-0#irrelevant since it's 1 stage
      trial_args$errtol<-0.01
      trial_args$maxpts<-10000
      trial_args$abseps<-.01

      full_case_args <- c(trial_args, case)
      
      if(trial_method=='cov') buildFun <- buildTrial
      if(trial_method=='MB') buildFun <- buildTrial_Maurer_Bretz_2013

      case$performance <- do.call(buildFun,full_case_args)$performance

      # check that power constraints are satisfied
      # if not, stop searching at this alpha configuration, and declare that we've failed to satisfy the requirements.
      # If no alpha configurations succeed, this will be the last time we change `all_constraints_satisfied`
      case_constraints_satisfied <- check_power_case(case)
      if(!case_constraints_satisfied) {
        all_constraints_satisfied <- FALSE
        break #no need to check the remaining cases
      }
    }

    #After we've searched over all stages, if we've met all constraints, stop the function.
    if(all_constraints_satisfied){
      return(c(
        'feasible'=1,
        'H01_eff_allocated'=alpha_1,
        'H02_eff_allocated'=alpha_2,
        'H0C_eff_allocated'=(FWER-alpha_1-alpha_2)
      ))
    }
    #print(trial_output$performance)
  }
}
#After we've searched over all configurations with no success, stop the trial.
return(c('feasible'=0))
}


# internal function
check_power_case <-function(case){

  power_ind <- c('Pow_H01','Pow_H02','Pow_H0C') #need to make sure we access powers in the same order for required and achieved power.
  power_mins <- unlist(case$power_mins[power_ind])
  power_achieved <- case$performance[power_ind]

  if(any(power_achieved < power_mins)) {
    constraints_satisfied <- FALSE
  }else{
    constraints_satisfied <- TRUE
  }

  return(constraints_satisfied)
}

#' binary search with arbitrary resolution
#'
#' The \code{\link{binsearch}} function in the \code{gtools} package searches only over integers. This is a wrapper for \code{\link{binsearch}} that also allows searching over the grid of points with distance \code{tol} between them. A target must also be entered (see \code{\link{binsearch}})
#' @param fun a function that determines the output over which we search (passed to \code{\link{binsearch}})
#' @param  tol resolution of the grid over which to search for \code{target} (see \code{\link{binsearch}})
#' @param  range a range over which to search for the input to \code{fun}
#' @param  ... passed to \code{\link{binsearch}}
#' @export
#' @import gtools
#' @seealso \code{\link{binsearch}}
#' @examples
#'  # best solution here is at x0,
#'  # which we can find with increasing precision
#'  x0 <- 10.241
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=1.00)
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=0.10)
#'  binsearchtol( function(x) x-x0, target=0, range=c(0,2*x0) , tol=0.05)
#'
binsearchtol <- function(fun, tol=1, range, ...){
  funTol <- function(x) fun(x*tol)

  soln <- binsearch(fun=funTol, range=range/tol, ...)
  soln$where <- soln$where*tol

  soln
}


#' Find the smallest total sample size such that a multi-stage trial meets desired constraints
#'
#' Takes as input a list of arguments (args) that define an adaptive trial (see \code{\link{buildTrial}} or \code{\link{optimizeTrial}}). This function adjusts the \code{n_total} argument in order to find the smallest maximum sample size that meets the power constraints specified in the \code{cases} argument.
#'
#' This function requires that the objective function contain a 'base' element, and a 'power_diffs' element that is nonnegative when power constraints are met. For example, see \code{\link{min_E_SS_power_constraints}}.
#'
#' @param args a list containing a subset of the arguments for the functions \code{\link{getEffBounds}} and \code{\link{simTrial}} (or comparable functions, see \code{trial_method} argument). This should include a FWER constraint.
#' @param cases A list of power constraints, in the same format as those sent to \code{\link{optimizeTrial}}
#' @param trial_method either 'cov' or 'MB' for Maurer-Bretz (2013).
#' @param objective_fun see \code{\link{optimizeTrial}}
#' @param min_n The smallest sample size to consider
#' @param max_n The largest sample size to consider
#' @param step_n The step size to consider when carrying out the binary search. For example, if \code{step_n = 10}, \code{min_n=0}, and \code{max_n=100}, the function will find the smallest n_total satisfying the supplied constraints, and that is also a multiple of 10.
#' @param showiter passed to \code{\link{binsearch}}
#' @export
#' @references
#' Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
#' @seealso \code{\link{min_n_feasible}}, \code{\link{feasibility_check}}
#' @return A list containing
#' \item{n}{The smallest feasible n_total}
#' \item{soln}{Output from \code{\link{binsearch}}}.
min_n_multistage<-function(args, cases, trial_method, objective_fun, min_n=1, max_n=min_n*1000, step_n=10, showiter=FALSE){


  fun4search<-function(x){
    args_x <- args
    args_x$n_total <- x
    args_x <- cleanArgList(args_x,max_K=args_x$n_total)

    perf_list<-get_case_perf_obj(base_args=args_x,
        cases = cases,
        skip_penalty=FALSE,
        objective_fun=objective_fun,
        trial_method=trial_method)
    # get_perf_mat(perf_list)

    # output here generally increasing after power constraint is met
    # and constant before that
    constraints_satisfied <- all(perf_list$obj$power_diffs>=0)
    if( constraints_satisfied) return(perf_list$obj$base)
    if(!constraints_satisfied) return(0)
    #if statements required, or we get Inf * 0.
  }

  soln<-binsearchtol(fun=fun4search,
    range=c(min_n,max_n),
    tol=step_n,
    target=.5,
    showiter=showiter
  )
  
  return(list(
    'n'=soln$where[length(soln$where)],
    'soln'=soln
    ))
}


#' @rdname feasibility_check
#' @export
#' @param min_n smallest sample size to consider for a 1-stage trial
#' @param max_n largest sample size to consider
#' @param step_n the step size for sample size. For example, if \code{step=20} then we will search in increments of 20 people.
#' @param ... passed to \code{feasibility_check}
min_n_feasible <- function(min_n=100, max_n = min_n * 10,step_n=5,showiter=FALSE, trial_method,...){
 
  # In order to use binsearch, we need to insert a monotonic function
  # We use the sample size multiplied by whether power constraints are met
  # This is zero for n_total too small, positive for n_total big enough, and monotonically increasing.
  # Our target is 0.5, which will highlight the point at which our
  # monotonic function just begins to be positive.

  last_fc <- 0 #used for it's names later on
  fc_tracker <- data.frame()
  fun4search <- function(x){   

    fc<-feasibility_check(n_total=x, trial_method=trial_method, ...)

    if(fc[1]==1){
      last_fc <<- fc
      fc_tracker <<- rbind(c('n_total'=x,fc),fc_tracker)
    }
    x*(fc[1]==1)
  }
  
  soln<-binsearchtol(fun=fun4search,
    range=c(min_n,max_n),
    tol=step_n,
    target=.5,
    showiter=showiter
    )
  names(fc_tracker)<-c('n_total',names(last_fc))

  out_n <- soln$value[2]

  if(all(soln$value==0)){
    warning('No solution found, max_n is too small. Returning max_n.')
    return(max_n)
  }
  if(soln$flag=='Lower Boundary'){
    warning('Solution found at min_n. Consider decreasing min_n, as a smaller trial might also be feasible.')
    out_n <- soln$value[1]
  }

  out <- (fc_tracker[fc_tracker$n_total==out_n,])[1,]
  return(unlist(out))
}


#' @rdname feasibility_check
#' @export
#' @param step_multiplier the step size
#' @param showiter passed to \code{\link[gtools]{binsearch}}, determines whether to show progress.
max_power_feasible <- function(n_total,cases, p1, trial_method,step_multiplier=.01,showiter=FALSE,FWER=0.025,npoints_sqrt=10,...){
 
  #The same approach is taken here as in min_n_feasible

  #check inputs
  if(step_multiplier<0 | step_multiplier>1) stop('step_multiplier must be between 0 and 1')


  #Check initial/default n_total performance
  fc1<-feasibility_check(n_total=n_total, p1, trial_method=trial_method,cases=cases, FWER=FWER,npoints_sqrt=npoints_sqrt)
  if(fc1[1]==1){
    warning('No search is necessary, initial n_total can achieve these power constraints. Returning original cases.')
    return(list(
      'multiplier'=1,cases=cases
      ))
  }

  #function to decrease min power required.
  deflate_mins<-function(z,multiplier){
    lapply(z,function(q) q*multiplier)
  }
  deflate_cases<-function(cases,multiplier){
    cases_2<-cases
    for(p in 1:length(cases)){
      cases_2[[p]]$power_mins <-
        deflate_mins(cases_2[[p]]$power_mins,multiplier=multiplier)
    }
    cases_2
  }

  #Search for largest multiplier meeting power constraints
  soln<-binsearchtol(fun=function(x){
     
      if(x==1){

        # recycle fc1 here to reduce redundant computation and potential
        # for inconsistencies due to random monte carlo error
        fc<-fc1

      }else{

        #multiply all power constraints by x
        cases_2<-deflate_cases(cases, multiplier=x)
        #check feasibility of resulting constraints
        fc<-feasibility_check(n_total=n_total, FWER=FWER, p1=p1, trial_method=trial_method,cases=cases_2, ...)

      }

      out <- ((1/(x+1))+10)* ( fc[1] == 1)
      out
      #divide by x+1 to avoid dividing by zero
      #For x <= optimal multiplier, out is decreasing (towards target=5) in x but bounded above 10
      #For x >  optimal multiplier, out is zero
      #We search for a target at 5
    },
    range=c(0,1),
    tol=step_multiplier,
    target=5.0,
    showiter=showiter
    )

  m <- soln$where[1]
  cases_final <- deflate_cases(cases,multiplier=m)
  return(list(
    multiplier = m,
    cases=cases_final
  ))
}


