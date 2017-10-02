# Functions to create and evaluate adaptive designs based on the covariance matrix

#Things to add:
# Decision theoretic problem - give reccomendations for who should get the treatment. Maximize # who benefit and minimize # who don't benefit.




#' Get smooth parametric boundaries for efficacy and futility
#'
#' A general parametric form for allocating alpha or creating futility boundaries. For details, see Fisher and Rosenblum (2016).
#' @param n_per_stage vector of sample sizes, or a vector proportional to these sample sizes
#' @param delta_vec a vector with elements greater than or equal to 0, one for each boundary to be computed (e.g. one for H_01, one for H_02, and one for H_0C).
#' @param coefficients_vec a vector of constraints, one for each boundary to be computed (e.g. one for H_01, one for H_02, and one for H_0C)
#' @param intercepts_vec a vector of constants to be added to the set to zero for efficacy boundaries
#' @param takeDiffs set to TRUE if calculating efficacy boundaries
#' @return A matrix with each row corresponding to one boundary (one hypothesis).
#' @references
#' Aaron Fisher and Michael Rosenblum (2016). Stochastic Optimization of Adaptive Enrichment Designs for Two Subpopulations. http://biostats.bepress.com/jhubiostat/paper279/
#' @export
#' @examples
#' getSmoothBounds(
#' 	  n_per_stage=1:5,
#' 	  intercepts_vec=c(0,0,1,1,0,1),
#' 	  delta_vec=c(1,1,1,1,1,1),
#' 	  coefficients_vec=c(1,1,1,1,2,2),
#'	  takeDiffs=FALSE
#' 	  )
getSmoothBounds <- function(n_per_stage, delta_vec, coefficients_vec, intercepts_vec, takeDiffs){
	
	K<-length(n_per_stage)
	lcv <- length(coefficients_vec)
	ldv <- length(delta_vec)
	liv <- length(intercepts_vec)
	H<-max(lcv,liv,ldv)

	if(any(delta_vec<0)) stop('delta must be nonnegative')
	
	if(any(c(lcv,liv,ldv)!=H)){
		warning('coefficients, delta, and intercepts are not the same length. All are being extended to match the maximum length.')
		if(lcv==1) coefficients_vec <- rep(coefficients_vec,H)
		if(ldv==1) delta_vec <- rep(delta_vec,H)
		if(liv==1) intercepts_vec <- rep(intercepts_vec,H)
	}

	out<-matrix(NA,H,K)
	
	if(takeDiffs & any(intercepts_vec!=0)) stop('Unexpected combination of takeDiffs and intercepts_vec')
	for(i in 1:H){
		x_base <- (cumsum(n_per_stage)/sum(n_per_stage))^delta_vec[i]
		increments <- x_base 
		if(takeDiffs) increments <- diff(c(0,x_base)) #used for efficacy boundaries as opposed to futility boundaries
		out[i,]<- intercepts_vec[i] + coefficients_vec[i] * increments 
	}

	return(out)
}




#x_null & x_len and is are vectors telling whether any of the arguments passed to a function are null, and how long they are.
checkEfficacyOverlap<-function(x_null, x_len){

	if(!x_null['delta_eff']){ #if delta_eff is entered
		if(any(!x_null[c( # checks whether any of these are specified
			'H01_eff_allocated',
			'H02_eff_allocated',
			'H0C_eff_allocated')])
			){

			stop("If delta_eff is entered, or vectors of efficacy allocated must not be entered")
		}
	}else{ #if delta_eff is *not* entered
		if(any(!x_null[c(# checks whether any of these are specified
		'H01_eff_total_allocated',
		'H02_eff_total_allocated',
		'H0C_eff_total_allocated')])
		){

			stop("If delta_eff is not entered, total efficacy allocated do not have interpretation.")
		}	
	}

}

checkFutilityOverlap<-function(x_null, x_len){
	if(!x_null['delta_futility']){ #if delta_futility is entered
		if(any(!x_null[c( # checks whether any of these are specified
			'H01_futility_boundaries',
			'H02_futility_boundaries',
			'H0C_futility_boundaries')])
			){

			stop("If delta_futility is entered, vectors of futility boundaries must not be entered.")
		}
	}else{ #if delta_futility is *not* entered
		if(any(!x_null[c(# checks whether any of these are specified
		'H01_futility_boundary_const',
		'H02_futility_boundary_const',
		'H0C_futility_boundary_const')])
		){

			stop("If delta_futility is not entered, futility coefficients_vec have no interpretation.")
		}	
	}


}


#' Generate efficacy boundaries and then calculate trial performance
#' This function first constructs the efficacy boundaries 
#' for an adaptive enrichment design by calling \code{\link{getEffBounds}}
#' and then simulates the trial design by calling \code{\link{simTrial}.
#' It ensures that efficacy boundaries are computed with the same arguments used to evaluate the trial's performance.
#'
#' Optionally, the user can specifically input \code{all_efficacy_boundaries} (or input \code{H01_efficacy_boundaries}, \code{H02_efficacy_boundaries}, and \code{H0C_efficacy_boundaries}), and \code{\link{getEffBounds}} will not be called. However, in such cases, it is simpler to just use the \code{\link{simTrial}} function directly. 
#'
#'
#'
#'
#' @param ... passed to \code{\link{getEffBounds}} and \code{\link{simTrial}}.
#' @export
#' @return the return value of \code{\link{simTrial}}
#'
buildTrial<-function(...){
	#############
	# Construct efficacy boundaries from alpha allocations
	# Get the performance of a given design
	#############

	#####################
	nmc <- names(match.call())

	#If we're missing efficacy boundaries, fill them in.
	if( !('all_efficacy_boundaries' %in% nmc) &
		!(all(c(
			'H01_efficacy_boundaries',
			'H02_efficacy_boundaries',
			'H0C_efficacy_boundaries'
			) %in% nmc))
	){
		
		all_efficacy_boundaries<-getEffBounds(...)

		#to avoid redundancy, assign other arguments to null 
		return(simTrial(
			all_efficacy_boundaries=all_efficacy_boundaries,
			'H01_efficacy_boundaries'=NULL,
			'H02_efficacy_boundaries'=NULL,
			'H0C_efficacy_boundaries'=NULL, 
			...
		))
	}
	#####################

	simTrial(...)

}



#' Compute efficacy stopping boundaries for an adaptive enrichment trial design based on
#' asymptotic, multivariate normal distribution (also called canonical distribution) of test statistics.
#' The result strongly controls the familywise Type I error rate, based on the 
#' generalized error-spending approach that allocates alpha (Type I error)
#' across stages and populations using the M_{COV} multiple testing procedure from:
#' Rosenblum, M., Qian, T., Du, Y., and Qiu, H., Fisher, A. (2016) Multiple Testing Procedures for Adaptive Enrichment Designs: Combining Group Sequential and Reallocation Approaches. Biostatistics. 17(4), 650-662. https://goo.gl/c8GlcH 
#' The algorithm for efficacy boundary construction involves sequential computation
#' of the multivariate normal distribution using the package mvtnorm.
#' Let \eqn{H01}, \eqn{H02} and \eqn{H0C} respectively denote the null hypotheses that there is no treatment effect in subpopulation 1, subpopulation 2 and the combined population.
#' 
#' @param p1 proportion of population in subpopulation 1.
#' @param r1 probability of being randomized to treatment in subpopulation 1
#' @param r2 probability of being randomized to treatment in subpopulation 2
#' @param var_s1_trt variance of the outcome under treament in subpopluation 1.
#' @param var_s1_con variance of the outcome under control in subpopluation 1.
#' @param var_s2_trt variance of the outcome under treament in subpopluation 2.
#' @param var_s2_con variance of the outcome under control in subpopluation 2.
#' @param time_limit time limit for calculations
#' @param num_stages number of stages for the trial
#' @param n_per_stage a vector with length equal to \code{num_stages}, telling the number of patient's outcomes to be observed in each stage. When there is no delay, this is equal to the number of patients enrolled per stage. When there is delay, this vector is not equal to the number of patients enrolled per stage.
#' @param n_total the total, maximum number of patients to recruit by the end of the study. If entered, n_per_stage will be scaled to have this sum.
#' @param FWER Familywise Type I error rate for the trial.
#' @param H01_eff_allocated a vector of length \code{num_stages} telling the proportion of Type I error to allocate to hypothesis \eqn{H01} at each stage of the trial.
#' @param H02_eff_allocated a vector of length \code{num_stages} telling the proportion of Type I error to allocate to hypothesis \eqn{H02} at each stage of the trial.
#' @param H0C_eff_allocated a vector of length \code{num_stages} telling the proportion of Type I error to allocate to hypothesis \eqn{H0C} at each stage of the trial.
#' @param FWER_allocation_matrix a matrix telling the proportion of Type I error to allocation to each hypothesis at each stage. If entered, this will override \code{H01_eff_allocated}, \code{H02_eff_allocated}, and \code{H0C_eff_allocated}.
#' @param H01_eff_total_allocated rather than setting the error allocated to each stage, the user can instead set the total error allocated to each hypothesis. \code{H01_eff_total_allocated}, \code{H02_eff_total_allocated}, and \code{H0C_eff_total_allocated} respectively tell the total Type I error to be allocated to \eqn{H01}, \eqn{H02}, and \eqn{H0C}. If set by the user, this will override the \code{H01_eff_allocated} vector.
#' @param H02_eff_total_allocated see \code{H01_eff_total_allocated}.
#' @param H0C_eff_total_allocated see \code{H01_eff_total_allocated}.
#' @param delta_eff This determines the allocation of Type I error across stages if \code{H01_eff_total_allocated}, \code{H02_eff_total_allocated} and \code{H0C_eff_total_allocated} are set by the user. See the source code.
#' @param abseps passed to pmvnorm in determining precision of calculations.
#' @param maxpts passed to pmvnorm in determining precision of calculations.
#' @param errtol determines precision of calculation of z-score boundary.
#' @param ... needed so that function ignores unused arguments when called by \code{\link{buildTrial}}
#' @export
#' @import mvtnorm 
#' @importFrom stats rnorm optim
#' @return a list of efficacy boundaries for the z-statistics corresponding to each null hypothesis.
#' @examples \dontrun{
#'
#' # Fully allocate the error for each stage
#' K <- 5
#' getEffBounds(p1 = 0.33,
#' 	 r1 = 1/2,
#' 	 r2 = 1/2,
#' 	 var_s1_trt = 0.375*(1-0.375),
#' 	 var_s1_con = 0.25*(1-0.25),
#' 	 var_s2_trt = 0.325*(1-0.325),
#' 	 var_s2_con = 0.2*(1-0.2),
#' 	 num_stages = 5,
#' 	 n_total = NULL,
#' 	 n_per_stage = rep(200,K),
#' 	 FWER = 0.025,
#'
#' 	 H01_eff_allocated=rep(0.025/(3*K),K),
#' 	 H02_eff_allocated=rep(0.025/(3*K),K),
#' 	 H0C_eff_allocated=rep(0.025/(3*K),K)
#' 	 )
#' # Use boundaries similar to O'Brien Flemming, with a parametric
#' # form that specifies the Type I error spent at each stage
#' # (note the changes in specifying the last four arguments)
#' getEffBounds(p1 = 0.33,
#' 	 r1 = 1/2,
#' 	 r2 = 1/2,
#' 	 var_s1_trt = 0.375*(1-0.375),
#' 	 var_s1_con = 0.25*(1-0.25),
#' 	 var_s2_trt = 0.325*(1-0.325),
#' 	 var_s2_con = 0.2*(1-0.2),
#' 	 num_stages = 5,
#' 	 n_total = NULL,
#' 	 n_per_stage = rep(200,K),
#' 	 FWER = 0.025,
#'
#' 	 delta_eff = .5, 
#' 	 H01_eff_total_allocated = 0.025/3,
#' 	 H02_eff_total_allocated = 0.025/3,
#' 	 H0C_eff_total_allocated = 0.025/3
#' 	)
#'
#'
#'}
getEffBounds<-function(p1,
	r1, #generally set to 1/2
	r2, #generally set to 1/2
	var_s1_trt,
	var_s1_con,
	var_s2_trt,
	var_s2_con,
	time_limit = 90,
	num_stages,
	n_total,
	n_per_stage,
	FWER, #= 0.025 generally
	H01_eff_allocated=NULL,
	H02_eff_allocated=NULL,
	H0C_eff_allocated=NULL,
	FWER_allocation_matrix=NULL,
	delta_eff=NULL,#set to 1 for approximately Pocock shaped boundaries
	H01_eff_total_allocated = NULL,
	H02_eff_total_allocated = NULL,
	H0C_eff_total_allocated = NULL,
	abseps,
	maxpts,
	errtol,
	...){

###
###
### Process for Computing Efficacy Boundaries
###
### First, construct cumulative sample size vectors
### Second, construct covariance matrix for statistics on z-scale
### Third, construct efficacy boundaries that correspond to alpha allocation
###
###

setTimeLimit(time_limit) # stops computation if taking greater than time_limit




if(!is.null(n_total)) n_per_stage <- n_total*n_per_stage/sum(n_per_stage)

p2 <- (1-p1)
outcome_variance_subpop_1 <- var_s1_trt/r1+var_s1_con/(1-r1)
outcome_variance_subpop_2 <- var_s2_trt/r2+var_s2_con/(1-r2)


ls_null<-sapply(ls(),function(x) is.null(eval(parse(text=x))))
ls_len<-sapply(ls(),function(x) length(eval(parse(text=x))))
checkEfficacyOverlap(x_null=ls_null, x_len=ls_len)

if(!is.null(delta_eff)){
	if(any(delta_eff < 0)){
		stop("Need nonnegative exponent")
	}

	eff_coeff_vec <- c(
	H01_eff_total_allocated,
	H02_eff_total_allocated,
	H0C_eff_total_allocated) #some of these may be NULL

	if(length(eff_coeff_vec) != 3 | any(is.na(eff_coeff_vec))){
		stop("If delta_eff is specified, all total efficacies allocated must also be specified.")
	}

	### Redo FWER allocation if using proportional to y^rho error spending function--we use delta_eff in place of rho; here y is the information accrued at a given analysis, which is proportional to the number of observed outcomes in our setup.

	FWER_allocation_matrix<-getSmoothBounds(n_per_stage=n_per_stage,
		coefficients_vec=eff_coeff_vec, delta_vec=delta_eff, intercepts_vec=rep(0,length(eff_coeff_vec)),takeDiffs=TRUE)
	rownames(FWER_allocation_matrix)<-c(
		'H01_eff_allocated',
		'H02_eff_allocated',
		'H0C_eff_allocated')

}else{
	FWER_allocation_matrix <- rbind(H01_eff_allocated,H02_eff_allocated,H0C_eff_allocated)
}



### Construct covariance matrix: subpopulation 1, subpopulation 2, combined population
covariance_matrix <- diag(3*num_stages)

# First set diagonal blocks corresponding to covariance of Z_{j,k} across stages for a given population j.
ss <- cumsum(n_per_stage) #Cumulative sample size in combined population
for(i in 1:num_stages){
for(j in 1:num_stages){
	covariance_matrix[i,j] <- 
	covariance_matrix[i+num_stages,j+num_stages] <-
	covariance_matrix[i+2*num_stages,j+2*num_stages] <- sqrt(min(ss[i],ss[j])/max(ss[i],ss[j]))
}}

# Next, set covariance of Z_{1,k},Z_{C,k'}
for(i in 1:num_stages){
for(j in 1:num_stages){
    covariance_matrix[i+2*num_stages,j] <-
    covariance_matrix[j,i+2*num_stages] <-
    sqrt((min(ss[i],ss[j])/max(ss[i],ss[j]))*(p1*outcome_variance_subpop_1/(p1*outcome_variance_subpop_1+p2*outcome_variance_subpop_2)))
 }}

# Last, set covariance of Z_{2,k},Z_{C,k'}
for(i in 1:num_stages){
for(j in 1:num_stages){
	covariance_matrix[i+2*num_stages,j+num_stages] <-
	covariance_matrix[j+num_stages,i+2*num_stages] <- sqrt((min(ss[i],ss[j])/max(ss[i],ss[j]))*(p2*outcome_variance_subpop_2/(p1*outcome_variance_subpop_1+p2*outcome_variance_subpop_2)))
 }}

### Construct efficacy boundaries from alpha allocations

ordering_of_statistics_by_stage <- as.vector(t(array(1:(3*num_stages),c(num_stages,3)))) ## Z_{1,1},Z_{2,1},Z_{C,1},...,Z_{1,K},Z_{2,K},Z_{C,K} where K=num_stages
covariance_matrix_ordered_by_stage <- covariance_matrix[ordering_of_statistics_by_stage,ordering_of_statistics_by_stage]
alpha_allocation <- as.vector(FWER*FWER_allocation_matrix/sum(FWER_allocation_matrix))

all_efficacy_boundaries <- c()
cumulative_alpha_allocation <- 0
for(index in 1:(3*num_stages)){
   cumulative_alpha_allocation <- cumulative_alpha_allocation + alpha_allocation[index]
   new_efficacy_boundary_upper_bound <- 20
   new_efficacy_boundary_lower_bound <- -20
   while(new_efficacy_boundary_upper_bound - new_efficacy_boundary_lower_bound > errtol){
    	new_efficacy_boundary_midpoint <- mean(c(new_efficacy_boundary_upper_bound,new_efficacy_boundary_lower_bound))
		cumulative_type_I_error <- 1-(pmvnorm(lower=rep(-Inf,index),upper=c(all_efficacy_boundaries,new_efficacy_boundary_midpoint),mean=rep(0,index),sigma=covariance_matrix_ordered_by_stage[1:index,1:index],algorithm=GenzBretz(abseps = abseps ,maxpts=maxpts))) 
		if(cumulative_type_I_error < cumulative_alpha_allocation){
			new_efficacy_boundary_upper_bound <- new_efficacy_boundary_midpoint
  		} else {new_efficacy_boundary_lower_bound <- new_efficacy_boundary_midpoint}
   }
   all_efficacy_boundaries <- c(all_efficacy_boundaries,new_efficacy_boundary_midpoint)
}

H01_efficacy_boundaries <- all_efficacy_boundaries[1+(3*(0:(num_stages-1)))]
H02_efficacy_boundaries <- all_efficacy_boundaries[2+(3*(0:(num_stages-1)))]
H0C_efficacy_boundaries <- all_efficacy_boundaries[3+(3*(0:(num_stages-1)))]

return(list(
	'H01_efficacy_boundaries'=H01_efficacy_boundaries,
	'H02_efficacy_boundaries'=H02_efficacy_boundaries,
	'H0C_efficacy_boundaries'=H0C_efficacy_boundaries
	))
}




#' Simulates an adaptive enrichment trial design to compute the following
#' performance criteria: power, expected sample size, and expected duration.
#' First, cumulative Z-statistics are constructed for each stage and population.
#' Next, the enrollment modification rule and multiple testing procedure are applied
#' at each stage, which determines when accrual is stopped for each subpopulation
#' and when (if at all) each population's null hypothesis is rejected. 
#' If efficacy boundaries have not yet been computed, the user should consider using \code{\link{buildTrial}}, which automatically completes this precursor step.
#'  
#' Let \eqn{H01}, \eqn{H02} and \eqn{H0C} respectively denote the null hypotheses that there is no treatment effect in subpopulation 1, subpopulation 2 and the combined population.
#'
#' @param p1 Proportion of population in subpopulation 1.
#' @param r1 probability of being randomized to treatment in subpopulation 1
#' @param r2 probability of being randomized to treatment in subpopulation 2
#' @param mean_s1_trt mean of the outcome under treament in subpopluation 1.
#' @param mean_s1_con mean of the outcome under control in subpopluation 1.
#' @param mean_s2_trt mean of the outcome under treament in subpopluation 2.
#' @param mean_s2_con mean of the outcome under control in subpopluation 2.
#' @param var_s1_trt variance of the outcome under treament in subpopluation 1.
#' @param var_s1_con variance of the outcome under control in subpopluation 1.
#' @param var_s2_trt variance of the outcome under treament in subpopluation 2.
#' @param var_s2_con variance of the outcome under control in subpopluation 2.
#' @param iter The number of simulated trials used to
#' estimate the power, expected sample size, and expected trial duration.
#' 
#' @param time_limit time limit for calculations.
#' @param n_per_stage a vector with length equal to \code{num_stages}, telling the number of patients to enroll in each stage.
#' @param n_total the total, maximum number of patients to recruit by the end of the study. If entered, n_per_stage will be scaled to have this sum.
#' @param num_stages
#' Total number of stages
#' used in each design  (\eqn{K}).  The maximum allowable number of stages is 20.
#' @param all_efficacy_boundaries a list of efficacy boundaries matching the output of \code{\link{getEffBounds}}
#' @param H01_efficacy_boundaries rather than setting \code{all_efficacy_boundaries}, the user can enter vectors for \code{H01_efficacy_boundaries}, \code{H02_efficacy_boundaries}, and \code{H0C_efficacy_boundaries}. 
#' @param H02_efficacy_boundaries see \code{H01_efficacy_boundaries}
#' @param H0C_efficacy_boundaries see \code{H01_efficacy_boundaries} 
#' @param H01_futility_boundaries a vector of futility boundaries for the hypothesis \eqn{H01}.
#' @param H02_futility_boundaries a vector of futility boundaries for the hypothesis \eqn{H02}.
#' @param H0C_futility_boundaries Not currently used in the algorithm, but may be added in the future.
#' @param delta_futility rather than setting the specific futility boundaries, parametric boundaries can be calculated. See \code{\link{getSmoothBounds}}.
#' @param H01_futility_boundary_const for use in \code{\link{getSmoothBounds}}
#' @param H02_futility_boundary_const for use in \code{\link{getSmoothBounds}}
#' @param H0C_futility_boundary_const for use in \code{\link{getSmoothBounds}}
#' @param enrollment_rate_combined The assumed
#' enrollment rate per year for the combined population.  This impacts the
#' expected duration of each trial design. Active enrollments from
#' the two subpopulations are assumed to be independent.  The enrollment rates
#' for subpopulations 1 and 2 are assumed proportional, based on \code{p_1}.
#' This implies that each stage of the adaptive design up to and including stage \code{k*} takes the same amount of time to complete, regardless of whether or not enrollment stops for subpopulation 2.  Each stage after \code{k*} will also take the same amount of time to complete. 
#' @param delay delay time from participant enrollment to observation	of his/her outcome (in years)
#' @param ... needed so that function ignores unused arguments when called by \code{\link{buildTrial}}
#' 
#' @details
#' This function is meant to be applied when there is prior
#' evidence that a treatment might work better in a one subpopulation
#' than in another. In this context, a trial with an adaptive enrollment
#' criteria would determine whether or not to continue enrolling patients
#' from each subpopulation based on interim analyses of whether each
#' subpopulation is benefiting. In order for the type I error and the
#' power of the trial to be calculable, the decision rules for changing
#' enrollment must be set before the trial starts. This function simulates trials with decision rules composed of efficacy boundaries and futility boundaries for \eqn{H01}, \eqn{H02}, and \eqn{H0C}, and reports the performance of the trial in
#' terms of power, expected sample size, and expected trial duration.
#'
#' @export
#'
simTrial <- function(
## Note: throughout, we denote the treatment arm by A=1 and control arm by A=0. 

## Subpopulation 1 proportion (Range: 0 to 1)
p1,
r1,
r2,

mean_s1_trt=NULL,
mean_s1_con=NULL,
mean_s2_trt=NULL,
mean_s2_con=NULL,

var_s1_trt =NULL,
var_s1_con =NULL,
var_s2_trt =NULL,
var_s2_con =NULL,

iter, 
time_limit = 90,
num_stages, 

n_total =  NULL,
n_per_stage,

all_efficacy_boundaries=NULL, #e.g. getEffBounds(). Arguments used here must match those used elsewhere in the function
H01_efficacy_boundaries=NULL, #null values of these vectors get filled in by the list.
H02_efficacy_boundaries=NULL,
H0C_efficacy_boundaries=NULL,

# Futility boundaries
H01_futility_boundaries=NULL, # Range (-10 to 10)
H02_futility_boundaries=NULL, # Range (-10 to 10)
H0C_futility_boundaries=NULL, # Range (-10 to 10)

#If set, these will override the above boundaries
delta_futility=NULL,
intercepts_futility=NULL,
H01_futility_boundary_const=NULL,
H02_futility_boundary_const=NULL,
H0C_futility_boundary_const=NULL,


# Enrollment rate for combined population (patients per year)
enrollment_rate_combined,
delay,
...
){




setTimeLimit(time_limit) # stops computation if taking greater than time_limit

if(!is.null(n_total)) n_per_stage <- n_total*n_per_stage/sum(n_per_stage)

p2 <- (1-p1)
outcome_variance_subpop_1 <- var_s1_trt/r1+var_s1_con/(1-r1)
outcome_variance_subpop_2 <- var_s2_trt/r2+var_s2_con/(1-r2)
SNR_subpop_1 <- (mean_s1_trt-mean_s1_con)/sqrt(outcome_variance_subpop_1)
SNR_subpop_2 <- (mean_s2_trt-mean_s2_con)/sqrt(outcome_variance_subpop_2)


#Replace any null efficacy boundaries with entries from the list
if(is.null(H01_efficacy_boundaries))
	H01_efficacy_boundaries<-all_efficacy_boundaries$H01_efficacy_boundaries
if(is.null(H02_efficacy_boundaries))
	H02_efficacy_boundaries<-all_efficacy_boundaries$H02_efficacy_boundaries
if(is.null(H0C_efficacy_boundaries))
	H0C_efficacy_boundaries<-all_efficacy_boundaries$H0C_efficacy_boundaries

## Override futility bounds with parametric bounds, if specified.


ls_null<-sapply(ls(),function(x) is.null(eval(parse(text=x))))
ls_len<-sapply(ls(),function(x) length(eval(parse(text=x))))
checkFutilityOverlap(x_null=ls_null, x_len=ls_len)

if(!is.null(delta_futility)){
	fut_coeff_vec <- c(
	H01_futility_boundary_const,
	H02_futility_boundary_const,
	H0C_futility_boundary_const) #some of these may be NULL

	if(length(fut_coeff_vec) != 3 | any(is.na(fut_coeff_vec)) ){
		stop("If delta_futility is specified, all futility boundary coefficients must also be")
	}
	if(length(intercepts_futility) != 3 | any(is.na(intercepts_futility)) ){
		stop("If delta_futility is specified, all futility boundary intercepts must also be")
	}

	fut_matrix <- getSmoothBounds(
		n_per_stage=n_per_stage,
		intercepts_vec=intercepts_futility,
		delta_vec=delta_futility,
		coefficients_vec=fut_coeff_vec,
		takeDiffs=FALSE
		)

	H01_futility_boundaries <- fut_matrix[1,]
	H02_futility_boundaries <- fut_matrix[2,]
	H0C_futility_boundaries <- fut_matrix[3,]

}



###
###
### Part II: Compute Design Performance
###
###

	cumulative_sample_size_vector_subpopulation_1 <- p1*cumsum(n_per_stage)
	cumulative_sample_size_vector_subpopulation_2 <- p2*cumsum(n_per_stage)
	# Enrollment rate subpop. 1 (patients per year)
	enrollment_rate_subpop_1 <- p1*enrollment_rate_combined
	# Enrollment rate subpop. 2 (patients per year)
	enrollment_rate_subpop_2 <- p2*enrollment_rate_combined 
        
## Get list of sample sizes corresponding to each interim analysis
	all_relevant_subpop_1_sample_sizes <- sort(unique(c(cumulative_sample_size_vector_subpopulation_1)))
	all_relevant_subpop_2_sample_sizes <- sort(unique(c(cumulative_sample_size_vector_subpopulation_2)))


## generate z-statistic increments (the change in the z-statistics at each stage)
	Z_subpop_1_increment <- array(0,c(length(all_relevant_subpop_1_sample_sizes),iter)) 
	Z_subpop_1_increment[1,] <- rnorm(iter)+SNR_subpop_1*sqrt(all_relevant_subpop_1_sample_sizes[1])
	if(length(all_relevant_subpop_1_sample_sizes)>1)
	{	for(i in 2:length(all_relevant_subpop_1_sample_sizes))
		{
			Z_subpop_1_increment[i,] <- rnorm(iter)+SNR_subpop_1*sqrt(all_relevant_subpop_1_sample_sizes[i]-all_relevant_subpop_1_sample_sizes[i-1])
		}
	}
	Z_subpop_2_increment <- array(0,c(length(all_relevant_subpop_2_sample_sizes),iter)) 
	Z_subpop_2_increment[1,] <- rnorm(iter)+SNR_subpop_2*sqrt(all_relevant_subpop_2_sample_sizes[1])
	if(length(all_relevant_subpop_2_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_2_sample_sizes))
		{
			Z_subpop_2_increment[i,] <- rnorm(iter)+SNR_subpop_2*sqrt(all_relevant_subpop_2_sample_sizes[i]-all_relevant_subpop_2_sample_sizes[i-1])
		}
	}
	
## generate partial sums of increments (weighted by per stage sample size)
## Construct cumulative z-statistics:
	# First for subpop_1 
	Z_subpop_1_partial_weighted_sum_of_increments <- Z_subpop_1_increment
	if(length(all_relevant_subpop_1_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_1_sample_sizes))
		{
			Z_subpop_1_partial_weighted_sum_of_increments[i,] <- 
		((sqrt(all_relevant_subpop_1_sample_sizes[i-1]/all_relevant_subpop_1_sample_sizes[i])*Z_subpop_1_partial_weighted_sum_of_increments[i-1,])		
			+ (sqrt((all_relevant_subpop_1_sample_sizes[i]-all_relevant_subpop_1_sample_sizes[i-1])/all_relevant_subpop_1_sample_sizes[i])*Z_subpop_1_increment[i,]))
		}
	}
	Z_subpop_1_cumulative <- array(0,c(num_stages,iter))
	for(i in 1:num_stages){
		index <- which(all_relevant_subpop_1_sample_sizes==cumulative_sample_size_vector_subpopulation_1[i])
		Z_subpop_1_cumulative[i,] <- Z_subpop_1_partial_weighted_sum_of_increments[index,]
	}
	# For subpopulation 2
	Z_subpop_2_partial_weighted_sum_of_increments <- Z_subpop_2_increment
	if(length(all_relevant_subpop_2_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_2_sample_sizes))
		{
			Z_subpop_2_partial_weighted_sum_of_increments[i,] <- 
		((sqrt(all_relevant_subpop_2_sample_sizes[i-1]/all_relevant_subpop_2_sample_sizes[i])*Z_subpop_2_partial_weighted_sum_of_increments[i-1,])		
			+ (sqrt((all_relevant_subpop_2_sample_sizes[i]-all_relevant_subpop_2_sample_sizes[i-1])/all_relevant_subpop_2_sample_sizes[i])*Z_subpop_2_increment[i,]))
		}
	}
	Z_subpop_2_cumulative <- array(0,c(num_stages,iter))
	for(i in 1:num_stages){
		index <- which(all_relevant_subpop_2_sample_sizes==cumulative_sample_size_vector_subpopulation_2[i])
		Z_subpop_2_cumulative[i,] <- Z_subpop_2_partial_weighted_sum_of_increments[index,]
	}
	# Define combined_population population z-statistics
	variance_component1 <- (p1^2)*outcome_variance_subpop_1/cumulative_sample_size_vector_subpopulation_1
	if(p2!=0){variance_component2 <- (p2^2)*outcome_variance_subpop_2/cumulative_sample_size_vector_subpopulation_2}else{variance_component2 <- 0*variance_component1}
	correlation_Z_subpop_1_with_Z_combined_population <- sqrt(variance_component1/(variance_component1+variance_component2)) 
	correlation_Z_subpop_2_with_Z_combined_population <- sqrt(variance_component2/(variance_component1+variance_component2))
	Z_combined_population_cumulative <- (correlation_Z_subpop_1_with_Z_combined_population*Z_subpop_1_cumulative + correlation_Z_subpop_2_with_Z_combined_population*Z_subpop_2_cumulative)
	
## Determine outcomes of each simulated trial

    # record if efficacy boundary ever crossed, for each of H0C and H01:
 	ever_cross_H0C_efficacy_boundary <- rep(0,iter)
	ever_cross_H01_efficacy_boundary <- rep(0,iter)
	ever_cross_H02_efficacy_boundary <- rep(0,iter)
	# indicator of stopping all enrollment, and of stopping only subpopulation 2, respectively:
    subpop_1_stopped <- rep(0,iter)
    subpop_2_stopped <- rep(0,iter)
    # indicators of rejecting null hypotheses:
	reject_H01 <- rep(0,iter)
	reject_H02 <- rep(0,iter)
	reject_H0C <- rep(0,iter)
    # record stage (just) after which enrollment stops for each subpopulation
	final_stage_subpop_1_enrolled_up_through <- rep(num_stages,iter)
        final_stage_subpop_2_enrolled_up_through <- rep(num_stages,iter)
	for(stage in 1:num_stages)
	{
          ever_cross_H0C_efficacy_boundary <- ifelse(Z_combined_population_cumulative[stage,]>H0C_efficacy_boundaries[stage],1,ever_cross_H0C_efficacy_boundary);
          ever_cross_H02_efficacy_boundary <- ifelse(Z_subpop_2_cumulative[stage,]>H02_efficacy_boundaries[stage],1,ever_cross_H02_efficacy_boundary)
          ever_cross_H01_efficacy_boundary <- ifelse(Z_subpop_1_cumulative[stage,]>H01_efficacy_boundaries[stage],1,ever_cross_H01_efficacy_boundary)
	# Determine if any new events where a null hypothesis is rejected for efficacy:
          reject_H01 <- ifelse((!subpop_1_stopped) & Z_subpop_1_cumulative[stage,]>H01_efficacy_boundaries[stage],1,reject_H01)
          reject_H02 <- ifelse((!subpop_2_stopped) & Z_subpop_2_cumulative[stage,]>H02_efficacy_boundaries[stage],1,reject_H02)
          reject_H0C <- ifelse((reject_H01 & reject_H02) | ((!subpop_1_stopped) & (!subpop_2_stopped) & Z_combined_population_cumulative[stage,]>H0C_efficacy_boundaries[stage]),1,reject_H0C)          
          subpop_1_stopped <- ifelse(reject_H01 | (Z_subpop_1_cumulative[stage,]<H01_futility_boundaries[stage]) | (Z_combined_population_cumulative[stage,]<H0C_futility_boundaries[stage]),1,subpop_1_stopped)    
          subpop_2_stopped <- ifelse(reject_H02 | (Z_subpop_2_cumulative[stage,]<H02_futility_boundaries[stage]) | (Z_combined_population_cumulative[stage,]<H0C_futility_boundaries[stage]),1,subpop_2_stopped)
       	 # record at what stage each subpop. stopped
         final_stage_subpop_1_enrolled_up_through <- ifelse((final_stage_subpop_1_enrolled_up_through==num_stages) & (subpop_1_stopped==1),stage,final_stage_subpop_1_enrolled_up_through)
         final_stage_subpop_2_enrolled_up_through <- ifelse((final_stage_subpop_2_enrolled_up_through==num_stages) & (subpop_2_stopped==1),stage,final_stage_subpop_2_enrolled_up_through)
	 }

pipeline_participant_max_subpopulation_1 <- p1*delay*enrollment_rate_combined
pipeline_participant_max_subpopulation_2 <- p2*delay*enrollment_rate_combined
max_sample_size_subpopulation_1 <- cumulative_sample_size_vector_subpopulation_1[num_stages]
max_sample_size_subpopulation_2 <- cumulative_sample_size_vector_subpopulation_2[num_stages]

#Distributions for sample size & duration
SS_dist<- pmin(
	cumulative_sample_size_vector_subpopulation_1[final_stage_subpop_1_enrolled_up_through] +
	pipeline_participant_max_subpopulation_1,max_sample_size_subpopulation_1)+
	pmin(cumulative_sample_size_vector_subpopulation_2[final_stage_subpop_2_enrolled_up_through]+
		pipeline_participant_max_subpopulation_2,
		max_sample_size_subpopulation_2)

dur_dist<- pmax(
	cumulative_sample_size_vector_subpopulation_1[final_stage_subpop_1_enrolled_up_through]/(ifelse(p1==0,Inf,p1)*enrollment_rate_combined),
	cumulative_sample_size_vector_subpopulation_2[final_stage_subpop_2_enrolled_up_through]/(ifelse(p2==0,Inf,p2)*enrollment_rate_combined)) +
	delay

return(list(performance=c(
'E_SS'=mean(SS_dist),
'E_dur'=mean(dur_dist), # expected duration
'Pow_H0C'=mean(reject_H0C), # power to reject H0C
'Pow_H01'=mean(reject_H01), # power to reject H01
'Pow_H02'=mean(reject_H02), # power to reject H02
'Pow_H01_and_H0C'=mean(reject_H0C & reject_H01), # power to reject H01 and H0C
'Pow_H02_and_H0C'=mean(reject_H0C & reject_H02), # power to reject H02 and H0C
'Pow_all'=mean(reject_H01 & reject_H02), # power to reject all (since H0C automatically rejected whenever both H01, H02 rejected)
'Pow_any'=mean(reject_H01 | reject_H02 | reject_H0C)), # power to reject at least one null hyp
'H01_efficacy_boundaries'=H01_efficacy_boundaries, # efficacy boundaries for null hypothesis H01
'H02_efficacy_boundaries'=H02_efficacy_boundaries, # efficacy boundaries for null hypothesis H02
'H0C_efficacy_boundaries'=H0C_efficacy_boundaries, # efficacy boundaries for null hypothesis H0C
'all_efficacy_boundaries'=all_efficacy_boundaries,
'H01_futility_boundaries'=H01_futility_boundaries, # futility boundaries for null hypothesis H01
'H02_futility_boundaries'=H02_futility_boundaries, # futility boundaries for null hypothesis H02
'H0C_futility_boundaries'=H0C_futility_boundaries, # futility boundaries for null hypothesis H0C
SS_dist = SS_dist, #Full distribution of SS
dur_dist = dur_dist #Full distribution of trial duration
))
	 
}




