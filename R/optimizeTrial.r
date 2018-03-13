# Functions for searching for an optimal design 
# in terms of E_SS or E_dur, subject
# to constraints on power & FWER

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".")) #ignore magrittr dot when doing R checks




#   ______                      __      __
#  /_  __/________ _____  _____/ /___ _/ /____
#   / / / ___/ __ `/ __ \/ ___/ / __ `/ __/ _ \
#  / / / /  / /_/ / / / (__  ) / /_/ / /_/  __/
# /_/ /_/   \__,_/_/ /_/____/_/\__,_/\__/\___/


# helper functions for translating between the 
# vector that optim needs and the list that users see.




# Convert vectorized arguments (from `unlist`) into the format expected by MR's functions for building & simulating trials.
# x is a vector of the total list of all possible arguments, after passing through `unlist`.
#Example for testing this function:
#x<-unlist(c(args_list_fixed,args_list_init))
#' @importFrom methods formalArgs
vector2list <- function(x){
  
  nx<-names(x)
  #argument names from MR's functions
  #note, check trial_method if making changes later on !!
  allArgNames <- unique(
    c(formalArgs(getEffBounds), 
    formalArgs(simTrial),
    formalArgs(getEffBounds_Maurer_Bretz_2013)) #include graph edges also
  )
  allArgNames<- allArgNames[allArgNames!='...']
  #All arguments in the list of arguments are either scalars or vectors.
  
  out<- list()
  for(i in 1:length(allArgNames)){
    arg_ind<-grepl(allArgNames[i],nx) #!! Future note - arg names shouldn't contain each other or grep will be confused
      ## !! FWER and FWER_allocation_matrix do, although it's not a problem now b/c we don't initilize and search for FWER_allocation_matrix
    if(!any(arg_ind)) next #skip un-specified args
    
    out[[allArgNames[i]]] <- x[arg_ind]
    
    #only keep high level names
    names(out[[allArgNames[i]]]) <- NULL
  }
  
  out
}

################ Clean arguments
#x is a *list* of arguments to all the functions in the buildTrial package
#We'll only be cleaning the list elements that are vectors, and ignoring the rest (which are redundant)
#test with:
# x<-c(args_list_fixed,args_list_init)
#Make sure to clean all of:
#unique(formalArgs(simTrial),formalArgs(getEffBounds))

###### helper functions 
edit_null<-function(a,min_a=-Inf,max_a=Inf,digits=NA){
  if( !is.null(a) ){
    a<-max(a,min_a)
    a<-min(a,max_a)
    if(!is.na(digits)) a <- round(a,digits=digits)
    return(a)
  }
  return(NULL)
}

sapply_null<-function(X,...){
  if(is.null(X[1])) return(NULL)
  sapply(X,...)
}
invLogit <- function(x){
  ex <- exp(x)
  out <- ex/(1+ex)

  out[which(is.infinite(ex) & ex > 0)] <- 1
  out
}
logit <- function(p){
  if(any( p>1 | p<0 )) stop('p must be between 0 and 1')
  log(p/(1-p))
}
######

#convert arguments to MR's functions to logit search space
#then, within the optimization objective function, convert back.
#x is a list of arguments
args_to_logit_space<-function(x,logit_search){
  out <- x #output will be the same except for at logit_search entries
  for(i in which(names(x)%in%logit_search)){
    out[[i]]<-logit(x[[i]])
  }
  out
}
logit_space_to_args<-function(x,logit_search){
  out<-x #output will be the same except for at logit_search entries
  for(i in which(names(x)%in%logit_search)){
    out[[i]]<-invLogit(x[[i]])
  }
  out
}

###### main cleaning function
#' Clean a list of arguments (for internal use)
#'
#' For internal use within scripts, and within \code{\link{optimizeTrial}}
#'
#' @export
#' @param x list of arguments
#' @param max_K maximum number of stages
#' @param min_n_total minimum per-stage sample size
cleanArgList<-function(x, max_K, min_n_total=10){ 
  
  smallNum <- 10^-10

  # (optional) Note for how this could be made cleaner:
    # Replace all instances of K or 1:K in this function with keep_ind. In addition to dropping stages >K, drop stages with n_stage=0.
      # This change would not quite be the same as our submitted work

  #First, we need to know how long our vectors should be.
  #We will use this throughout the function.
  K<- edit_null(x$num_stages, min_a=1, max_a=max_K,digits=0)
  if(is.null(x$num_stages)) K<-max_K #max_K shouldn't ever be null, since we require it as part of "core_args"
  x$num_stages <- K
  
  
  x$n_per_stage <- sapply_null(x$n_per_stage, function(z){
    edit_null(z,min_a=0) #this is later standardized to have sum n_total
  })[1:K] 
  x$n_per_stage[1] <- edit_null(x$n_per_stage[1],min_a=smallNum)
    #zero people in stage 1 is both uninterpretable and will cause NaNs in covariance matrix later on 

  
  x$n_total <- edit_null(x$n_total,min_a=min_n_total)
  if(!is.null(x$n_total) & !is.null(x$n_per_stage)){
    x$n_per_stage <-  x$n_total * x$n_per_stage[1:K] / sum(x$n_per_stage[1:K])
  }
  x$enrollment_rate_combined <-
    edit_null(x$enrollment_rate_combined, min_a=5)
  
  ##############
  # Code in this chunk can also be handled by searching
  # over the logit space for nonnegative or decimal parameters.
  # So, if also using logit_search, some of these checks may be redundant.
  
  x$p1 <- edit_null(x$p1,min_a=smallNum, max_a = 1-smallNum)
  x$r1 <- edit_null(x$r1,min_a=smallNum, max_a = 1-smallNum)
  x$r2 <- edit_null(x$r2,min_a=smallNum, max_a = 1-smallNum)
  x$FWER <- edit_null(x$FWER,min_a=smallNum, max_a = 1-smallNum)
  
  if(!is.null (x$delta_eff)) x$delta_eff <- sapply(x$delta_eff, function(z) edit_null(z,min_a=0))
  if(!is.null (x$delta_futility)) x$delta_futility <- sapply(x$delta_fut, function(z) edit_null(z,min_a=0))

  x$H01_eff_allocated <- sapply_null(x$H01_eff_allocated, 
    function(z) edit_null(z,min_a=smallNum))[1:K]
  x$H02_eff_allocated <- sapply_null(x$H02_eff_allocated,
    function(z) edit_null(z,min_a=smallNum))[1:K]
  x$H0C_eff_allocated <- sapply_null(x$H0C_eff_allocated,
    function(z) edit_null(z,min_a=smallNum))[1:K]
  
  x$H01_eff_total_allocated <- edit_null(x$H01_eff_total_allocated, min_a=smallNum) #!! smallNum is used to avoid NaNs in get Eff bounds
  x$H02_eff_total_allocated <- edit_null(x$H02_eff_total_allocated, min_a=smallNum)
  x$H0C_eff_total_allocated <- edit_null(x$H0C_eff_total_allocated, min_a=smallNum)

  if( (!is.null(x$H01_eff_total_allocated)) &
      (!is.null(x$H02_eff_total_allocated)) &
      (!is.null(x$H0C_eff_total_allocated))){
         sum_eff_totals<- x$H01_eff_total_allocated + x$H02_eff_total_allocated + x$H0C_eff_total_allocated
         x$H01_eff_total_allocated <- x$H01_eff_total_allocated/sum_eff_totals
         x$H02_eff_total_allocated <- x$H02_eff_total_allocated/sum_eff_totals
         x$H0C_eff_total_allocated <- x$H0C_eff_total_allocated/sum_eff_totals
  }
   

  x$graph_edge_12 <- edit_null(x$graph_edge_12, min_a=0, max_a=1)
  x$graph_edge_2C <- edit_null(x$graph_edge_2C, min_a=0, max_a=1)
  x$graph_edge_C1 <- edit_null(x$graph_edge_C1, min_a=0, max_a=1)
  ##############
  
  
  
  x$H01_futility_boundaries <- x$H01_futility_boundaries[1:K]
  x$H02_futility_boundaries <- x$H02_futility_boundaries[1:K]
  x$H0C_futility_boundaries <- x$H0C_futility_boundaries[1:K]
  
  #sapply_null is used for var_**_*** when a discrete prior is entered.
  if(!is.null(x$var_s1_trt)) x$var_s1_trt <- sapply_null(x$var_s1_trt,function(z)
    max(z,.000000001)
  )
  if(!is.null(x$var_s1_con)) x$var_s1_con <- sapply_null(x$var_s1_con,function(z)
    max(z,.000000001)
  )
  if(!is.null(x$var_s2_trt)) x$var_s2_trt <- sapply_null(x$var_s2_trt,function(z)
    max(z,.000000001)
  )
  if(!is.null(x$var_s2_con)) x$var_s2_con <- sapply_null(x$var_s2_con,function(z)
    max(z,.000000001)
  )
  
  #######
  # Adjust allocated errors so they sum to FWER
  # Note, this must come *after* we truncate efficacy allocated by K.
  # This is already done within the code for getEffBounds, but
  # nice to have also in the clean output from optimizeTrial
  if(is.null(x$FWER_allocation_matrix) & 
      !any(is.null(x$H01_eff_allocated)) &
      !any(is.null(x$H02_eff_allocated)) &
      !any(is.null(x$H0C_eff_allocated)) )
    x$FWER_allocation_matrix <- as.matrix(rbind(x$H01_eff_allocated,x$H02_eff_allocated,x$H0C_eff_allocated)[,1:K])
  if(!is.null(x$FWER_allocation_matrix) & !is.null(x$FWER)){ #a second is.null essentially checks if all of the information to compute this was included in x
    x$FWER_allocation_matrix <- as.matrix(x$FWER_allocation_matrix[,1:K])
    x$FWER_allocation_matrix <- x$FWER * 
      as.matrix(x$FWER_allocation_matrix[,1:K])/
      sum(x$FWER_allocation_matrix[,1:K])
    x$H01_eff_allocated <- x$FWER_allocation_matrix[1,1:K]
    x$H02_eff_allocated <- x$FWER_allocation_matrix[2,1:K]
    x$H0C_eff_allocated <- x$FWER_allocation_matrix[3,1:K]
  }
  #######	
  
  
  #######
  #It is generally not reccomended to vary these with optim
  x$iter<- edit_null(x$iter,min_a=100,digits=0)
  x$time_limit <- edit_null(x$time_limit, min_a=60)
  #######
  
  x
}
################


















#     ______            __            __
#    / ____/   ______ _/ /_  ______ _/ /____
#   / __/ | | / / __ `/ / / / / __ `/ __/ _ \
#  / /___ | |/ / /_/ / / /_/ / /_/ / /_/  __/
# /_____/ |___/\__,_/_/\__,_/\__,_/\__/\___/

# Functions for evaluating a design, 
# and calculating the objective function


#' Returns the expected sample size (mostly for internal use)
#'
#' Mostly for internal use. This function is called by \code{\link{get_case_perf_obj}}.
#' @param cases Describes the scenarios under which sample size should be calculated. Each case should include a \code{performance} element and a \code{weight} element.
#' @export
#' @return the expected sample size across cases.
get_E_SS<-function(cases){
  sum(unlist(lapply(cases,
                    function(x){
                      x$performance['E_SS'] * x$weight}
  )))
}

# skip_penalty: sometimes valid designs will incur a penalty due to monte 
# carlo error. If you know, in theory, that a design shouldn't have a
# penalty, you can push it not to with this argument. This is equivalent 
# to using base_obj, but stores the result in a different 
# way (for useful for consistency with other version of this code).
min_perf_power_constraints<-function(cases, metric, skip_penalty, exp_p=3){
  ########
  
  #Extract power from each case we care about
  achieved_power <- 
  power_mins <-  matrix(NA,length(cases),3) 
  
  testchar<-c('Pow_H01','Pow_H02','Pow_H0C') #use index vector to ensure the same order between required power and achieved power
  for(p in 1:length(cases)){
    achieved_power[p,] <- cases[[p]]$performance[testchar]
    power_mins[p,] <- unlist(cases[[p]]$power_mins[testchar])
  }
  if(any(c(power_mins)>=1)) stop('Minimum power must be between 0 and 1.')
  
  #######
  #Calculate objective function to minimize:
  power_penalty <- (abs(power_mins-achieved_power)*100)^exp_p * (achieved_power < power_mins) 
    #Let the penalty be pth order differentiable
  if(skip_penalty) power_penalty <- 0 
    #!! In future, could adjust this and functions that 
    # access history of objective function to just use base_obj, and 
    # not bother with this skip_penalty
  base_obj <- sum(unlist(lapply(cases, function(x){
      x$performance[metric] * x$weight
  })))
  penalized_obj <- sum(power_penalty) + base_obj
  
  ######
  names(penalized_obj)<-c()

  # additional diagnostics to assess the found solution  
    # negative indicates insufficient power.
  power_diffs <- achieved_power - power_mins

  return(list(penalized = penalized_obj, base = base_obj, power_diffs = power_diffs))
}


#' Sample objective function to minimize
#'
#' \code{\link{optimizeTrial}} requires users to input an objective function to evaluate the trial design in each case. \code{min_E_SS_power_constraints} is a pre-set option for an objective function based on expected sample size. \code{min_E_dur_power_constraints} is a pre-set option for an objective function based on expected duration. Both functions' output a list including a penalized objective function value. Any custom objective function supplied to  \code{\link{optimizeTrial}} must be able to run in the context of the \code{\link{get_case_perf_obj}} function, and give output in the same format as that of the functions below. 
#' @param ... for internal use.
#' @export
#' @return a list containing
#' \item{base}{the expected sample size, or expected duration}
#' \item{penalized}{the value of \code{base} plus a penalty if any of the power constraints are violated}
#' \item{power_diffs}{a matrix showing the difference between achieved power and required power for each case, and each constraint.}
min_E_SS_power_constraints <- function(...) {
  min_perf_power_constraints(metric = 'E_SS',...)
}


#' @rdname min_E_SS_power_constraints
#' @export
min_E_dur_power_constraints <- function(...) {
  min_perf_power_constraints(metric = 'E_dur',...)
}


#' Get trial performance for each case. 
#' 
#' This function is mostly for internal use within \code{\link{optimizeTrial}}. It first gets trial performance for each case. Once the characteristics of the trial in each case are simulated, the objective function can be calculated.
#'
#' See argument definitions in \code{\link{optimizeTrial}} for more information.
#'
#' @param ... is passed to \code{objective_fun}
#' @param base_args the trial design and other inputs to functions determined by \code{trial_method}. This corresponds to the combined list of \code{args_list_init} and \code{args_list_fixed} from \code{\link{optimizeTrial}}.
#' @param cases analogous to the input for \code{\link{optimizeTrial}}.
#' @param objective_fun a function to evaluate the design, see \code{\link{optimizeTrial}}.
#' @param trial_method either 'MB' or 'cov', depending on the type of design being used. See \code{\link{optimizeTrial}}
#' @param return_entries a list of strings telling which elements of the returned value from \code{objective_fun} should be stored. Set to 'all' to return all elements from functions determined by \code{trial_method}.
#' 
#' @return A list with elements
#' \item{cases}{a modified version of \code{cases} in which each element has been appended with descriptive metrics from functions determined by \code{trial_method}.}
#' \item{obj}{The output from \code{objective_fun}}
#'
#' @export
get_case_perf_obj<-function(base_args, cases, objective_fun, trial_method, return_entries=c('performance'),...){
  
  
  #get output from simulating power & SS cases
    for(i in 1:length(cases)){

      all_args_i <- c(cases[[i]],base_args)

      #!! Note for future improvements: In all of these cases,
      #  is it possible to instead get the formalArgs
      #  of all sub-functions called using `...`?
      if(trial_method=='cov'){
        arg_names_i <-intersect(
          names(all_args_i),
          c(formalArgs(getEffBounds), formalArgs(simTrial))
        )
        buildTrialMethod <- buildTrial
      }

      if(trial_method=='MB'){
        arg_names_i <- intersect(
          names(all_args_i),
          c(formalArgs(getEffBounds_Maurer_Bretz_2013),
            formalArgs(simTrial_Maurer_Bretz_2013))
        )
        buildTrialMethod <- buildTrial_Maurer_Bretz_2013
      }

      builtTrial <- do.call(buildTrialMethod, all_args_i[arg_names_i])
      for(j in 1:length(builtTrial)){
        nj <- names(builtTrial)[j]
        if(nj %in% return_entries | all(return_entries=='all')){
          cases[[i]][nj]<-builtTrial[nj]
        }
      }
      

    }
  
  ####### Objective function, and saving results
  objs <- objective_fun(cases = cases, ...)
  
  return(list('obj'=objs, cases=cases))
}


#' Get matrix with trial design performance in each supplied case
#' 
#' This function parses the output of \code{\link{get_case_perf_obj}} into a more easily read form.
#' @param perf_list output object from \code{\link{get_case_perf_obj}}
#' @export
get_perf_mat<-function(perf_list){
  t(as.data.frame(
    lapply(perf_list$cases,function(x){
      x$performance
    })
  ))
}


















###################################
###################################
###################################
###################################
#                __  _           _
#   ____  ____  / /_(_)___ ___  (_)___  ___
#  / __ \/ __ \/ __/ / __ `__ \/ /_  / / _ \
# / /_/ / /_/ / /_/ / / / / / / / / /_/  __/
# \____/ .___/\__/_/_/ /_/ /_/_/ /___/\___/
#     /_/
#
# User facing function to search for an approximately optimal design


#' Use simulated annealing to search over a space of trial designs 
#'
#' We aim to minimize an objective function which we aim to minimize is the expected sample size, with additive penalties on power. If power for any alternative hypothesis is less than a supplied threshold (see \code{cases}), a penalty will be applied.
#'
#' The argument \code{trial_method} determines which type of trial should be performed ('cov' or 'MB'). Based on this, the relevant trial parameters to be optimized should be listed in \code{args_list_init}. The parameters to hold fixed should be listed in \code{args_list_fixed}.
#'
#' The \code{cases} argument should be a list of cases in which power, expected sample size and/or expected duration should be calculated. Internally, at each iteration of the search, the \code{\link{get_case_perf_obj}} function is used to evaluate the trial.
#'
#' @param args_list_init a list containing a subset of the arguments for the functions \code{\link{getEffBounds}} and \code{\link{simTrial}} (or comparable functions, see \code{trial_method argument}). This is the subset of arguments over which we will search when trying to minimize the objective function. The values supplied in \code{args_list_init} form the initial values at which we start the search. Elements of this list must be vectors.
#' @param args_list_fixed a subset of arguments for the functions \code{\link{getEffBounds}} and \code{\link{simTrial}} (or comparable functions, see \code{trial_method argument}) which we will not search over. Instead, these arguments will remain fixed at the user-specified values. Arguments not specified will be set to their defaults. Elements of this list must be vectors.
#' @param cases a list of lists describing power constraints, and prior distribution of treatment effects over which to calculate the expected sample size (or duration). List elements are individually passed to \code{objective_fun} for evaluation (via \code{\link{get_case_perf_obj}}).
#' @param max_K an upper limit for the number of stages for the trial. \code{optimizeTrial} will not search for trials with more stages than this.
#' @param objective_fun objective function to minimize. This depends on the evaluation of the current proposed design at each element of \code{cases}. Pre-set options that can be used are \code{\link{min_E_SS_power_constraints}} and \code{\link{min_E_dur_power_constraints}}. If a custom objective function is supplied here instead, it must be able to run in the context of the \code{\link{get_case_perf_obj}} function, and give output in the same format as that of \code{\link{min_E_SS_power_constraints}}.
#' @param logit_search a list of parameters for which optim will search in the logit space. This can include parameters which are  bounded to the range (0,1).
#' @param local_n_search [DEPRECATED] (logical) whether a search over n_total should be done at each iteration of optim. This had previously generated the output \code{local_search_calls}.
#' @param max_n_local [DEPRECATED] largest sample size for searches within iterations of optimization. Especially useful to set for optimizing a 1-stage design.
#' @param min_n_local [DEPRECATED] smallest sample size for searches within iterations of optimization. Especially useful to set for optimizing a 1-stage design.
#' @param time_limit_optim_cpu used to set a session time limit (see \code{\link{setSessionTimeLimit}}).
#' @param time_limit_optim_elapsed used to set a session time limit (see \code{\link{setSessionTimeLimit}}).
#' @param max_n_feas largest total sample size to consider when checking feasibility
#' @param optim_method passed to \code{\link{optim}}.
#' @param maxit passed to \code{\link{optim}}. Tells the number of optimization iterations. For simulated annealing, this is the same as the number of function calls. 
#' @param trial_method the type of trial to run. 'cov' for covariance-based, 'MB' for Maurer Bretz (2013).
#' @param parscale_ratio_N used to create the \code{parscale} argument passed to \code{\link{optim}}. Tells the change in sample size that is comparable to a unit change in other parameters.
#' @param parscale_ratio_K used to create the \code{parscale} argument passed to \code{\link{optim}}. Tells the change in the number of stages that is comparable to a unit change in other parameters.
#' @param parscale_ratio_logit used to create the \code{parscale} argument passed to \code{\link{optim}}. For any parameter where do search for optimal values on the logit space, this tells the unit change that is comparable to a unit change in other parameters.
#' @param build_precision if TRUE, the search will start by making approximate calculations for the efficacy boundaries, and will increase the precision of these efficacy boundary calculations as the search progresses.
#' @param print_interval if set to \code{NULL}, no intermediate results will be reported. If set to a positive integer, this integer tells the number of simulated annealing steps to run in between results being printed to the console. The last number printed in each line of results is best objective function yet found. This report is a more detailed version of the results generally reported by \code{\link{optim}}, if \code{control$trace=1}.
#' @param npoints_sqrt passed to \code{\link{feasibility_check}} to see if it is indeed possibly to have a trial at this sample size that meets the power constraints. If not, no optimization is performed.
#' @param stage1_feasible output from \code{\link{min_n_feasible}}, if pre-calculated. If supplied, trial feasibility will be assumed.
#' @param returnFullPath (logical) tells whether the entire history of the search should be returned as well.
#' @param verbose (logical) should progress be printed via \code{\link{message}}.
#' @param print_obj (logical) In progress reports, should the objective function be printed rather than just the interpretable, unpenalized expected sample size.
#' @return A list of results with
#'
#' \item{soln}{the best design found during optimization}
#' \item{performance_crossvalidated}{The solution found by optim may have low sample size or duration only due to Monte-Carlo error in calculating trial performance. To counter this potential bias, we recalculate the performance of the design returned by our optimization, and return the results in \code{performance_crossvalidated}.}
#' \item{obj_final_crossvalidated}{The objective function value of the design from optim. Like \code{performance_crossvalidated}, this value is recalculated after the optimization procedure.}
#' \item{optim_results}{output from \code{\link{optim}}}
#' \item{obj_penalized_trajectory}{Each evaluation of the objective function.}
#' \item{fullPath}{design and results from each iteration of the optimization}
#' \item{time_optimized}{time taken for optimization procedure}
#' \item{stage1_feasible}{output from \code{\link{min_n_feasible}}}
#' \item{feasible}{[DEPRECATED]}
#' \item{optim_iterations}{Number of iterations performed}
#' \item{local_search_calls}{[DEPRECATED]}
#'
#' @export
#' @references Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
#' @import mvtnorm dplyr GenSA
#' @examples \dontrun{
#' 
#' 
#' 
#' set.seed(0)
#' 
#' ## Generate example inputs
#' inputsMISTIE <- getExampleInputsMISTIE(iter=50) # ~ 30 seconds
#' str(inputsMISTIE,1) # list of inputs
#' 
#' 
#' 
#' ## Optimize adaptive trial
#' optimized <- do.call(optimizeTrial,inputsMISTIE) #~ 4-7 minutes
#' 
#' 
#' ## Visualize results for approximately optimized trial
#' soln <- optimized$soln
#' 
#' 
#' eff_bounds<-getEffBoundFromOptimSoln(
#'   soln = soln,
#'   case = inputsMISTIE$cases[[1]]
#'   )
#' fut_bounds <- cbind(
#'     soln$H01_futility_boundaries,
#'     soln$H02_futility_boundaries,
#'     soln$H0C_futility_boundaries
#'   )[(1:soln$num_stages-1),] 
#' 
#' matplot(data.frame(eff_bounds),type='o',pch=1:3,lty=1:3,
#'   ylim=range(c(eff_bounds,fut_bounds)),
#'   ylab='Boundary (z-scale)',
#'   xlab='Cumulative number with outcome observed',
#'   col='blue', main='Trial decision boundaries')
#' matlines(fut_bounds, col='red',type='o',pch=1:3,lty=1:3)
#' legend('bottomright',
#'   c('H01 Eff','H02 Eff', 'H0C Eff','H01 Fut','H02 Fut', 'H0C Fut'),
#'   col=rep(c('blue','red'),each=3),
#'   pch=rep(1:3,times=2),lty=rep(1:3,times=2))
#' 
#' }
optimizeTrial<-function(
  args_list_init, 
  args_list_fixed, 
  cases,
  max_K = 20, 
  objective_fun = min_E_SS_power_constraints,  
  logit_search = c(),
  local_n_search=FALSE, 
  max_n_local=NULL,
  min_n_local=NULL,
  time_limit_optim_cpu=NA,
  time_limit_optim_elapsed=NA,
  max_n_feas=NULL,
  optim_method = 'SANN',
  maxit = 1000, 
  trial_method,
  parscale_ratio_N = 100, 
  parscale_ratio_K = 1, 
  parscale_ratio_logit = 3, 
  build_precision = TRUE, 
  print_interval=NULL,
  npoints_sqrt=25,
  stage1_feasible=NULL,
  returnFullPath=FALSE,
  verbose = is.null(print_interval),
  print_obj = FALSE
){
  
  
  
  
  ###########
  # Check certain inputs

  # Check that core arguments are inputted
  core_args<-c('num_stages','n_per_stage','n_total','FWER','p1','r1')
  args_list_all<-c(args_list_init,args_list_fixed)
  missing_core_args<-core_args[!core_args %in% names(args_list_all)]
  if(length(missing_core_args)>0) stop(paste(c('Must include the following variables in either args_list_fixed or args_list_init:',missing_core_args),collapse=' '))
  if(any(duplicated(names(args_list_all)))) stop('args_list_init and args_list_fixed cannot share element names')
  if(args_list_all['p1']>=1 | args_list_all['p1'] <=0) stop('p1 must be set betwen 0 and 1') # Possibly redundant with other checks, but that is good. E_dur will be set to NaN in some cases if this happens.
  if('FWER_allocation_matrix' %in% args_list_all) stop('optimization of FWER_allocation_matrix not supported. Use individual alpha allocations instead.') #to avoid the issue with FWER and FWER_allocation_matrix being contained in each other and confusing grep.

  #Omit parts of logit search that we don't search over, as this can
  #inverse transform the wrong arguments later on
  logit_search <- intersect(logit_search, names(args_list_init))
  

  #Check to make sure weights of expected sample size cases sum to 1.
  case_weights <- unlist(lapply(cases, function(x){
    if(is.null(x$weight)) stop('Weights (possibly 0) must be defined for cases')
    if(is.na(x$weight)) stop('Weights (possibly 0) must be defined for cases')
    x$weight
  }), use.names=FALSE)
  if(any(case_weights < 0)) stop('Weights for cases must be nonnegative')
  if(sum(case_weights) != 1){
    warning('Standardizing weights for cases to sum to 1')
    case_weights<- case_weights/ sum(case_weights)
    for(i in 1:length(cases)){
      cases[[i]]$weight <- case_weights[i]
    }
  }
    
  #strip away extra element names from list inputs, as these will mess up indexing later on.
  for(i in 1:length(args_list_init)){
    names(args_list_init[[i]]) <- NULL
  }
  for(i in 1:length(args_list_fixed)){
    names(args_list_fixed[[i]]) <- NULL
  }

  ###########
  
  ###########
  #Check trial feasibility 
  #Code here answers: What will be the first value we use for n_total? (on iteration 1 of the search?)
  n_total_iter1 <- args_list_init$n_total
  if(is.null(n_total_iter1)) n_total_iter1 <- args_list_fixed$n_total
  
  search4n<- 'n_total'%in%names(args_list_init)
  
  if(verbose & is.null(stage1_feasible))
    message('Checking feasibility of power constraints.\nAs a reference, searching for the smallest 1-stage trial that meets the same power constraints...')
  
  
  n_feasible<-Inf #We always need to calculate this, as something to compare our final (optimized) design to.
  n_search_warn<- tryCatch({

    if(is.null(max_n_feas)) max_n_feas <- max(n_total_iter1,1) * 100 #!!hard-coded!!
    stage1_args_list <- list()
    if(trial_method!='cov') stage1_args_list <- list( graph_edge_12=1/2,graph_edge_2C=1/2, graph_edge_C1=1/2)
    if(is.null(stage1_feasible)){
      stage1_feasible <- min_n_feasible(
        min_n = 10,
        max_n = max_n_feas,
        step_n=1,
        trial_method=trial_method,
        p1=args_list_all$p1,
        r1=args_list_all$r1,
        showiter=FALSE,
        FWER=args_list_fixed$FWER,
        cases=cases,
        npoints_sqrt=npoints_sqrt,
        trial_args = stage1_args_list
      )
    }
    n_feasible<-stage1_feasible['n_total']
  
  },warning=function(x) x)
  
  
  if(n_total_iter1 < n_feasible){
    warning_part1<-paste0('Even with a 1-stage trial, it does not appear to be possible to meet these power constraints at this',c(' initial')[search4n],' maximum sample size (n_total).\nTo increase precision of feasibility check, increase npoints_sqrt.\n')
    
    if(verbose) message(c('Initial ')[search4n],'n_total not feasible, searching for alternate n_total value and/or alternate power constraints...')
    
    ######
    if('warning' %in% class(n_search_warn)){
      n_feasible <- NA
      warning_part2<-'No feasible sample size found'
      power_feasible = list(cases=NA, multiplier=NA)
    }else{
      if(search4n){
        n_total_iter1 <-
          args_list_init$n_total <- n_feasible*1.25
        message('Minimum required n_total was found at ',
                n_feasible,'. Resetting initial n_total to ',args_list_init$n_total,
                ' (1.25 factor increase above minimum).')
        warning_part2<-''
        
      }
      if(!search4n){
        power_feasible <-max_power_feasible(
          n_total=n_total_iter1,
          trial_method=trial_method,
          FWER=args_list_fixed$FWER,
          cases=cases,
          npoints_sqrt=npoints_sqrt,
          step_multiplier=.01,
          showiter=FALSE,
          p1=args_list_all$p1,
          r1=args_list_all$r1
        )
        
        warning_part2<-paste0('To achieve a feasible design, consider setting n_total at least as high as ',
                              n_feasible,', or reducing power constraints (cases) by a factor of ',
                              power_feasible$multiplier)

      }
    } 
    ######
    
    warning(paste(warning_part1, warning_part2))
    
    if( (!search4n) | 'warning' %in% class(n_search_warn) ){
      return(list('feasible'=FALSE,
                  'n_feasible'=n_feasible,
                  'power_feasible' = power_feasible
      ))
    }
    
  }
  
  if(verbose) message('Feasibility checked. Beginning search for optimal trial...\n')
  
  ###########
  
  
 
  #starting value of theta
  #First convert the appropriate parameters to the logit space.
  theta_init<-args_to_logit_space(args_list_init,logit_search) %>%
    unlist
  
  
  
  
  
  # optim_iterations, fullPath and 
  # local_search_calls (local_search_calls is now DEPRECATED) 
  # are defined in this environment, later on, 
  # just before calling optim
  
  # Function (of theta) which optim will call and optimize.
  # This function calls back up to it's enclosing environment to get args_list_fixed.
  fun2opt<-function(theta, return_full_output=FALSE){

    #Need to add this line if using GenSA, since GenSA strips the names off of theta before using it.
    names(theta)<-names(theta_init) 

    #Track how many times the function has been called (!!redundant with optim counts output).
    optim_iterations <<- optim_iterations+1

    ########
    # Get list of arguments for trial building & simulating functions
    vec_fixed <- args_to_logit_space(args_list_fixed,logit_search) %>% unlist
    if(any(duplicated(c(names(vec_fixed),names(theta))))){
      stop('duplicate arguments across fixed and variable parameters')
    }
    
    args_all <- c(vec_fixed,theta) %>% 
      vector2list %>% 
      logit_space_to_args(.,logit_search) %>%
      cleanArgList(., max_K=max_K, min_n_total=n_feasible)

    # Adjust arguments that vary within the optim procedure
    if(build_precision){
      args_all$abseps <-  0.001 / min(1,(1000*optim_iterations/maxit))
      args_all$maxpts <- 1000 * max(1,(100*optim_iterations/maxit))
      #Changing args_all$iter doesn't help speed very much
    }
    
    if(!local_n_search) local_search_calls <<- c(local_search_calls,0)
    if( local_n_search){
      # This option is no longer being used or maintained.
      # It is listed as DEPRECATED in the documentation
      # It did not appear to help optimization very much,
      # although it might yet be improved in future versions.
    
      
      if(is.null(max_n_local)) max_n_local <- max(n_feasible,args_all$n_total)*5      
      if(is.null(min_n_local)) min_n_local <- n_feasible
      ##!! Future note - hard-coding could be removed
      # User specified max and min allows this to be 
      # used with either a 1-stage or multi-stage optimization.
      
      suppressWarnings({
        n_local <-
        min_n_multistage(args_all,
          cases,
          trial_method,
          objective_fun,
          min_n=min_n_local,
          max_n=max_n_local,
          step_n=10,
          showiter=FALSE
        )
      })
        # warnings are possible if either boundary is hit, but sensible results will still be returned. (the upper boundary is always returned)
      
      local_search_calls <<- c(local_search_calls, n_local$soln$numiter)

      args_all$n_total <- n_local$n

    }

    #Evaluate the performance of the design for each case
    iter_time<-system.time({
      perf_list<-get_case_perf_obj(args_all,
                                cases = cases,
                                skip_penalty=FALSE,
                                objective_fun=objective_fun,
                                trial_method=trial_method) 
    })
    

    ##############################
    ##############################
    ##############################
    # Code chunk for storing intermediate values of objective function

    perf_mat<-get_perf_mat(perf_list)

    power_satisfied <- unlist(lapply(perf_list$cases, function(x) check_power_case(x)
    ))
    E_SS <- get_E_SS(perf_list$cases)

    updateFullPath(
      design=args_all,
      obj=perf_list$obj,
      E_SS=E_SS,
      cases=perf_list$cases,
      time= iter_time,
      power_satisfied=power_satisfied
    )
    
    #print intermediate results
    #!! Currently, this code chunk only works for finding "best" trials, not "worst" trials.
    if(!is.null(print_interval) & verbose){
      if(optim_iterations %% print_interval ==0 ){
        best_obj <- lapply(fullPath, function(z){z$obj$penalized}) %>%
          unlist %>% min(.,na.rm=TRUE)

        print_penalized_text<-NULL
        if(print_obj) print_penalized_text <- paste0(';  penalized obj=',round(perf_list$obj$penalized,2),' / ',round(best_obj,2))
        
        message(paste0('iter-',optim_iterations,': constraints_met = ',paste(power_satisfied,collapse=', '),print_penalized_text))
      }}
    
    ##############################
    ##############################
    ##############################
    
    if(return_full_output){
      perf_mat
      return(list('performance'=perf_mat,'obj'=perf_list$obj))
    }else{
      return(perf_list$obj$penalized)
    }
  }
  
  

  
  
  #Set up parscale (argument for optim)
  parscale <- rep(1,length(theta_init))
  names(parscale) <- names(theta_init)
  
  
  parscale[	c(grep('n_per_stage',names(parscale)), 
    grep('n_total',names(parscale)))  ] <- 
    parscale_ratio_N
  if('num_stages' %in% names(theta_init)){
    parscale['num_stages'] <- parscale_ratio_K
  }
  
  parscale[names(parscale) %in% logit_search] <- parscale_ratio_logit
  
  
  #################################
  # Optim implementation
  
  
  ######## 
  # Variables to keep track of the progress of optim.
  # These variables are defined in the parent env of fun2opt, but not the global env.
  fullPath<- list()
  optim_iterations <- 0 # A counter to see how many times fun2opt is called.
  local_search_calls <- c() #vector of iterations for local search at each iteration of optimization.
  
  #This function must be defined in same env as fullPath.
  updateFullPath<-function(design, obj, cases, time=0,
  E_SS= get_E_SS(cases),
  power_satisfied= NULL ){

    # Could store results more efficiently, but this is not the bottleneck in the code.
    # Generally using lists here rather than data-frames,
    # to allow possibility that theta is nonconstant length


    if(is.null(power_satisfied)){
      power_satisfied <- unlist(lapply(cases, function(x) check_power_case(x)
      ))
    }

    case_performance<-lapply(cases, function(x) x$performance)

    fullPath[[length(fullPath)+1]] <<- list(
      'design'=design,
      'obj'=obj,
      'E_SS'=E_SS,
      'case_performance'=case_performance,
      'time'= time,
      'power_satisfied'=power_satisfied
    )

  }


  #Set up the 1-stage feasible design to be the first element of fullPath
  n_feasible_unnamed<-n_feasible
  names(n_feasible_unnamed) <- c()
  vector_1stage<-c(
      'num_stages'=1,
      'n_per_stage'=n_feasible_unnamed,
      'H01_futility_boundaries'=-Inf,
      'H02_futility_boundaries'=-Inf,
      'H0C_futility_boundaries'=-Inf,
      'graph_edge_12'=1/2, #!! Need to specify, otherwise we get missing args for MB trials. It would be better to have as defined from a separate search, but this is just a rough starting point.
      'graph_edge_2C'=1/2,
      'graph_edge_C1'=1/2,
      stage1_feasible[names(stage1_feasible)!='feasible']
    )

  stage1_design <- c(
    unlist(args_list_fixed[!names(args_list_fixed)%in%names(vector_1stage)]),vector_1stage) %>% #no need to transform from logit space
    vector2list %>% 
    cleanArgList(.,max_K=1, min_n_total=n_feasible)

  perf_list_1stage <- get_case_perf_obj(stage1_design,
      cases = cases,
      skip_penalty = TRUE,
      objective_fun=objective_fun,
      trial_method=trial_method) 
  
  updateFullPath(
    design=stage1_design,
    obj=perf_list_1stage$obj,
    cases=perf_list_1stage$cases,
    time= 0,
    power_satisfied = rep(TRUE,length(cases)) #This might not be technically true due to rounding error, but for now we assume there is no rounding error.
  )
  ######## 

  control<-list(parscale=parscale)
  control$maxit<-maxit
  
  
  lower <- rep(-Inf,length(theta_init))
  upper <- rep(Inf,length(theta_init))
  
  #get names of arguments searched for on real line
  logit_search_names_vec <- names(unlist(args_list_all[logit_search]))
  nti <- names(theta_init)
  nti[nti %in% logit_search_names_vec] <- ''

  alpha_inds<- grepl('eff_allocated',nti) 
  lower[alpha_inds]<-0
  upper[alpha_inds]<-1
  lower[grepl('n_',nti)] <- 1
  lower[nti=='num_stages'] <- 2

  if(optim_method=='GenSA'){
    controlGenSA<-list(maxit=maxit)
  }

  optimized<-NULL 
  limit_time_cpu <- !is.na(time_limit_optim_cpu)
  limit_time_elapsed <- !is.na(time_limit_optim_elapsed)
  time_optimized<-system.time(try({
    if(limit_time_cpu)     setSessionTimeLimit(cpu = time_limit_optim_cpu) #reset to Inf after optim procedure
    if(limit_time_elapsed) setSessionTimeLimit(elapsed = time_limit_optim_elapsed) #reset to Inf after optim procedure  
    
    if(optim_method=='sann_random'){
      optimized<-sann_random(par=theta_init, fn=fun2opt, control=control)
    }else if(optim_method=='GenSA'){
      optimized<-GenSA(theta_init, fun2opt,lower=lower,
              upper=upper,control=controlGenSA)
    }else if(optim_method %in% c('L-BFGS-B','Brent')){
      optimized<-optim(theta_init, fun2opt, lower=lower,
              upper=upper, method=optim_method,control=control)
    }else{
      optimized<-optim(theta_init, fun2opt, method=optim_method,control=control)
    }
    
  }))
  if(limit_time_cpu) setSessionTimeLimit(cpu = Inf)
  if(limit_time_elapsed) setSessionTimeLimit(elapsed = Inf)
  if(limit_time_cpu | limit_time_elapsed){
    warning('When time limits are used (i.e. when either time_limit_optim_cpu or time_limit_optim_elapsed is not NA), the time limit is reset to Inf after the search calculations have completed, via setSessionTimeLimit')
  }

  #################################
  # Check & report results

  # Note - this will find the best solution even if optim timed out.
  previous_penalized_obj <- unlist(lapply(fullPath, function(z){z$obj$penalized}))

  solution_ind <- which(previous_penalized_obj==min(previous_penalized_obj))[1]
  soln <- fullPath[[solution_ind]]$design
  
  if(	(!is.null(args_list_fixed$FWER)) &
      solution_ind==1
  ){
    warning('No sequential trial was found with a better performance than that of a standard 1-stage trial testing the same hypothesis. Returning a trial with only 1 stage.')
  }else{
    # If we have found an admissable multistage trial design:
    # update n_total, such that trial has minimum n needed to achieve power
    # do this *before* cross-validated assessment of trial performance
    # otherwise, we run the risk of getting a trial with too low 
    # power, due to a chance occurence that an underpowered trial
    # happened to do well in the optimation seach iteration.

    # !! Note !! This step makes the procedure slightly less flexible
    # as it is assumed that power is monotonically
    # increasing in n_total (without penalties).

    n_total_crossval <-
      min_n_multistage(soln, cases, trial_method, objective_fun, min_n=soln$n_total/5, max_n=soln$n_total*5, step_n=1, showiter=FALSE)$n #!! hardcoded range

    soln$n_total <- n_total_crossval
    soln <- cleanArgList(soln, max_K = soln$num_stages)
  }
  
  #Cross-validation
  perf_list<-get_case_perf_obj(soln,
      cases = cases,
      objective_fun=objective_fun,
      skip_penalty=FALSE,
      trial_method=trial_method) 	
  
  perf_mat<-get_perf_mat(perf_list)
  
  if(!returnFullPath) fullPath <-NULL
  
  
  return(list(
    'soln'=soln, #We report the full solution here. Reporting just part of this list of arguments (e.g. just the solutions from optim) is tricky because in order to clean and process them we need all arguments anyway (e.g. we need FWER)
    'performance_crossvalidated' = perf_mat,
    'obj_final_crossvalidated' = perf_list$obj,
    'optim_results'= optimized,
    'obj_penalized_trajectory'=previous_penalized_obj, #Note, this is not strictly decreasing.
    'fullPath' = fullPath,
    'time_optimized'=time_optimized,
    'stage1_feasible'=stage1_feasible,
    'feasible'=TRUE,#!!deprecated!! (noted in documentation)
    'optim_iterations'=optim_iterations,
    'local_search_calls'=local_search_calls #!!deprecated!! (noted in documentation)
  ))
}


#' Get efficacy boundary from a solution found by \code{\link{optimizeTrial}}
#'
#' This function helps to parse the output from \code{\link{optimizeTrial}} into a more easily viewed format.
#' @export
#' @param soln the \code{soln} element in the returned list from \code{\link{optimizeTrial}}
#' @param case a scenario in which to calculate boundaries, in the same format as an element of the \code{cases} list argument for \code{\link{optimizeTrial}}
getEffBoundFromOptimSoln <- function(soln, case){
  args_all <- c(soln,case)
  ind_for_eff <- intersect(names(args_all),formalArgs(getEffBounds))
  args_eff <- args_all[ind_for_eff]

  do.call(getEffBounds,args_eff)
}

