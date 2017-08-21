# Functions to create and evaluate adaptive designs based on Maurer Bretz & alpha reallocation


#Things to add:
# Decision theoretic problem - give reccomendations for who should get the treatment. Maximize # who benefit and minimize # who don't benefit.

#' Generate efficacy boundaries and then calculate trial performance using method from Maurer Bretz (2013)
#'
#' This function calls \code{\link{getEffBounds_Maurer_Bretz_2013}} and \code{\link{simTrial_Maurer_Bretz_2013}}.
#'
#' @param ... passed to \code{\link{getEffBounds_Maurer_Bretz_2013}} and \code{\link{simTrial_Maurer_Bretz_2013}}.
#' @export
#' @references Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
buildTrial_Maurer_Bretz_2013 <-function(...){
	simTrial_Maurer_Bretz_2013(...)
}


#' Compute efficacy boundaries for an adaptive trial using method from Maurer Bretz (2013)
#'
#' All arguments are analagous to \code{\link{getEffBounds}}, with the addition of three parameters: \code{graph_edge_12}, \code{graph_edge_2C}, and \code{graph_edge_C1}.
#' 
#' @param graph_edge_12 The proportion of alpha to reallocate from H_01 to H_02 in the event that H_01 is rejected
#' @param graph_edge_2C The proportion of alpha to reallocate from H_02 to H_0C in the event that H_02 is rejected
#' @param graph_edge_C1 The proportion of alpha to reallocate from H_0C to H_01 in the event that H_0C is rejected
#'
#' @export
#' @import mvtnorm
#' @importFrom stats rnorm
#' @references Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
getEffBounds_Maurer_Bretz_2013 <-function(p1,
	r1,
	r2,
	var_s1_trt,
	var_s1_con,
	var_s2_trt,
	var_s2_con,
	time_limit,
	num_stages,
	n_total = NULL,
	n_per_stage,
	FWER,
	H01_eff_allocated=NULL,
	H02_eff_allocated=NULL,
	H0C_eff_allocated=NULL,
	FWER_allocation_matrix=NULL,
	delta_eff=NULL,
	H01_eff_total_allocated = NULL,
	H02_eff_total_allocated = NULL,
	H0C_eff_total_allocated = NULL,
	abseps,
	maxpts,
	errtol,
	graph_edge_12,
	graph_edge_2C,
	graph_edge_C1,
	...){



setTimeLimit(time_limit) # stops computation if taking greater than time_limit


ls_null<-sapply(ls(),function(x) is.null(eval(parse(text=x))))
ls_len<-sapply(ls(),function(x) length(eval(parse(text=x))))
checkEfficacyOverlap(x_null=ls_null, x_len=ls_len)

if(!is.null(n_total)) n_per_stage <- n_total*n_per_stage/sum(n_per_stage)

p2 <- (1-p1)
outcome_variance_subpop_1 <- var_s1_trt/r1+var_s1_con/(1-r1)
outcome_variance_subpop_2 <- var_s2_trt/r2+var_s2_con/(1-r2)

### Redo FWER allocation if using proportional to y^rho error spending function--we use delta_eff in place of rho; here y is the information accrued at a given analysis, which is proportional to the number of observed outcomes in our setup.

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

	FWER_allocation_matrix<-getSmoothBounds(n_per_stage=n_per_stage,
		coefficients_vec=eff_coeff_vec, delta_vec=delta_eff, intercepts_vec=rep(0,length(eff_coeff_vec)), takeDiffs=TRUE)
	rownames(FWER_allocation_matrix)<-c(
		'H01_eff_allocated',
		'H02_eff_allocated',
		'H0C_eff_allocated')


	H01_eff_allocated <- FWER_allocation_matrix[1,]
	H02_eff_allocated <- FWER_allocation_matrix[2,]
	H0C_eff_allocated <- FWER_allocation_matrix[3,]
	
}
FWER_allocation_matrix <- rbind(H01_eff_allocated,H02_eff_allocated,H0C_eff_allocated)

### Construct covariance matrix for z-statistics of a single population over time (it is same for all populations)
covariance_matrix <- diag(num_stages)
ss <- cumsum(n_per_stage) #Cumulative sample size in combined population
for(i in 1:num_stages){
for(j in 1:num_stages){
	covariance_matrix[i,j] <- sqrt(min(ss[i],ss[j])/max(ss[i],ss[j]))
}}

##
## Precompute boundaries for each possible subset of rejected null hypotheses under Maurer Bretz algorithm (Section 3.2 of their paper).
##

# Below, ordering of null hypothesis is H01, H02, H0C being first, second, third (only relevant in understanding the notation)
# Inputs to this algorithm are graph edges and weights (where each node in the graph is an elementary null hypothesis):
   	graph_edge_initial_values <- t(array(c(NA,graph_edge_12,1-graph_edge_12,1-graph_edge_2C,NA,graph_edge_2C,graph_edge_C1,1-graph_edge_C1,NA),c(3,3))) # (i,j) entry represents directed edge value from hypothesis i to hypothesis j, using above ordering
	hypothesis_initial_weights <- c(sum(H01_eff_allocated),sum(H02_eff_allocated),sum(H0C_eff_allocated))/sum(c(H01_eff_allocated,H02_eff_allocated,H0C_eff_allocated)) # used in Bonferroni adjustments

	all_efficacy_boundaries <- array(NA,c(2,2,2,3,num_stages)) # first three indices correspond to status (0=unrejected; 1=rejected) of each null hypothesis; fourth index represents null hypothesis/population under consideration; fifth represents stage. Each entry in the arry is the corresponding boundary for null hyp. under consideration after rejection of indicated subset of nulls at given stage. E.g., boundaries[0,0,1,2,] is the list of efficacy boundaries for subpopulation 2 (H02) after third null hyp (H0C) has been rejected.
#Search over all possible subsets of elementary null hypotheses (where 0 indicates not yet rejected; 1 indicates rejected):
	for(h1 in 0:1){	for(h2 in 0:1){	for(h3 in 0:1){ 
	       rejection_status_vector <- c(h1,h2,h3)
	       # initialize graph
	       graph_edge_values <- graph_edge_initial_values
	       hypothesis_weights <- hypothesis_initial_weights
	       # compute updated graph after rejections
	       for(null_hyp in 1:3){   
	       		if(rejection_status_vector[null_hyp]==1){
			     # if null hypothesis rejected, update graph accordingly using Algorithm from Maurer and Bretz, Section 3.2
				for(hyp_temp in 1:3){if(rejection_status_vector[hyp_temp]==0) {hypothesis_weights[hyp_temp] <- hypothesis_weights[hyp_temp] + graph_edge_values[null_hyp,hyp_temp]*hypothesis_weights[null_hyp]}}
				hypothesis_weights[null_hyp] <- 0;
				for(hyp_temp_1 in 1:3){for(hyp_temp_2 in 1:3){if(rejection_status_vector[hyp_temp_1]==0 && rejection_status_vector[hyp_temp_2]==0 && hyp_temp_1 != hyp_temp_2 && graph_edge_values[hyp_temp_1,null_hyp]*graph_edge_values[null_hyp,hyp_temp_1]<1){graph_edge_values[hyp_temp_1,hyp_temp_2] <- ((graph_edge_values[hyp_temp_1,hyp_temp_2] + graph_edge_values[hyp_temp_1,null_hyp]*graph_edge_values[null_hyp,hyp_temp_2])/(1-graph_edge_values[hyp_temp_1,null_hyp]*graph_edge_values[null_hyp,hyp_temp_1]))}}}
		 	}                    
	       }
	       #print(cbind(h1,h2,h3,graph_edge_values,hypothesis_weights))
	       # determine boundaries for each null hypothesis over all stages
	       for(null_hyp in 1:3){
	       		    alpha_allocation <- FWER*hypothesis_weights[null_hyp]*switch(null_hyp,H01_eff_allocated/sum(H01_eff_allocated),H02_eff_allocated/sum(H02_eff_allocated),H0C_eff_allocated/sum(H0C_eff_allocated))			   
	       		    if(is.na(sum(alpha_allocation))) stop('invalid alpha allocation')#AF - Note, I think this error check is because we could end up with 0/0 in switch, above? (!!)
	       		    efficacy_boundary_existing <- c()
			    cumulative_alpha_allocation <- 0
			    for(index in 1:num_stages){
			        cumulative_alpha_allocation <- cumulative_alpha_allocation + alpha_allocation[index]
			        new_efficacy_boundary_upper_bound <- 20
	   		        new_efficacy_boundary_lower_bound <- -20
   			        while(new_efficacy_boundary_upper_bound - new_efficacy_boundary_lower_bound > errtol){
    			           new_efficacy_boundary_midpoint <- mean(c(new_efficacy_boundary_upper_bound,new_efficacy_boundary_lower_bound))
			           cumulative_type_I_error <- 1-(pmvnorm(lower=rep(-Inf,index),upper=c(efficacy_boundary_existing,new_efficacy_boundary_midpoint),mean=rep(0,index),sigma=covariance_matrix[1:index,1:index],algorithm=GenzBretz(abseps = abseps ,maxpts=maxpts))) 
				   if(cumulative_type_I_error < cumulative_alpha_allocation){
				      new_efficacy_boundary_upper_bound <- new_efficacy_boundary_midpoint
  				   } else {new_efficacy_boundary_lower_bound <- new_efficacy_boundary_midpoint}
   			        }
				efficacy_boundary_existing <- c(efficacy_boundary_existing,new_efficacy_boundary_midpoint)
			    }
	       	       	    all_efficacy_boundaries[h1+1,h2+1,h3+1,null_hyp,] <- efficacy_boundary_existing 
		}
	}}}
	       #print("efficacy boundaries")
	       #for(h1 in 0:1){ for(h2 in 0:1){ for(h3 in 0:1){print(all_efficacy_boundaries[h1+1,h2+1,h3+1,,])}}}
	       #print("weights")
	       #print(hypothesis_initial_weights)
	return(all_efficacy_boundaries)
}

#' Simulate a trial to compute its power, expected sample size, and expected duration using method from Maurer Bretz (2013)
#'
#' All arguments are analagous to \code{\link{simTrial}}. However, additional parameters are required by \code{\link{getEffBounds_Maurer_Bretz_2013}}.
#'
#' @param ... passed to \code{\link{getEffBounds_Maurer_Bretz_2013}}
#'
#' @export
#' @references Maurer, W. and F. Bretz (2013). Multiple testing in group sequential trials using graphical approaches. \emph{Statistics in Biopharmaceutical Research.}
simTrial_Maurer_Bretz_2013 <- function(
## Note: throughout, we denote the treatment arm by A=1 and control arm by A=0. 

## Subpopulation 1 proportion (Range: 0 to 1)
p1,
r1,
r2,

mean_s1_trt,
mean_s1_con,
mean_s2_trt,
mean_s2_con,

var_s1_trt, 
var_s1_con, 
var_s2_trt, 
var_s2_con, 

iter, 
time_limit,
num_stages, 

n_total =  NULL,
n_per_stage,

all_efficacy_boundaries=NULL, 

# Futility boundaries
H01_futility_boundaries, # Range (-10 to 10)
H02_futility_boundaries, # Range (-10 to 10)
H0C_futility_boundaries, # Range (-10 to 10)

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

if(is.null(all_efficacy_boundaries)){	
    all_efficacy_boundaries <- getEffBounds_Maurer_Bretz_2013(p1=p1,r1=r1,r2=r2,var_s1_trt=var_s1_trt,var_s1_con=var_s1_con,var_s2_trt=var_s2_trt,var_s2_con=var_s2_con,num_stages=num_stages,n_per_stage=n_per_stage, time_limit=time_limit,...)
}

if(!is.null(n_total)) n_per_stage <- n_total*n_per_stage/sum(n_per_stage)

p2 <- (1-p1)
outcome_variance_subpop_1 <- var_s1_trt/r1+var_s1_con/(1-r1)
outcome_variance_subpop_2 <- var_s2_trt/r2+var_s2_con/(1-r2)
SNR_subpop_1 <- (mean_s1_trt-mean_s1_con)/sqrt(outcome_variance_subpop_1)
SNR_subpop_2 <- (mean_s2_trt-mean_s2_con)/sqrt(outcome_variance_subpop_2)



## Override futility bounds with parametric bounds, if specified.
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
	

##
## Determine outcomes of each simulated trial   
##

   
    # indicators of rejecting null hypotheses:
	reject_H01 <- rep(0,iter)
	reject_H02 <- rep(0,iter)
	reject_H0C <- rep(0,iter)
    # record stage (just) after which enrollment stops for each subpopulation
	final_stage_subpop_1_enrolled_up_through <- rep(num_stages,iter)
        final_stage_subpop_2_enrolled_up_through <- rep(num_stages,iter)

   for(simulated_trial_iteration in 1:iter){         
      ## For each simulated trial, keep track of if/when hypotheses rejected, and weights/edges of re-allocation graph
        stage_when_hypotheses_rejected <- c(0,0,0) # 0 indicates not yet rejected; ordering is H01, H02, H0C
	rejection_status_vector <- c(0,0,0)
	enrollment_status <- c(1,1,1) #indicators of subpop. 1 enrollment, subpop. 2 enrollment, and both enrollment up through current stage. 
   	for(stage in 1:num_stages){
	    new_null_hypothesis_rejected_current_stage <- TRUE
	    while(new_null_hypothesis_rejected_current_stage){
		new_null_hypothesis_rejected_current_stage <- FALSE
		for(null_hyp in 1:3){
		     if(rejection_status_vector[null_hyp]==0 && enrollment_status[null_hyp]==1){ # if null_hyp not yet rejected and enrollment not previously stopped, check if it can be rejected at current stage
			current_z_statistic <- switch(null_hyp,Z_subpop_1_cumulative[stage,simulated_trial_iteration],Z_subpop_2_cumulative[stage,simulated_trial_iteration],Z_combined_population_cumulative[stage,simulated_trial_iteration])		   
			if(current_z_statistic > all_efficacy_boundaries[rejection_status_vector[1]+1,rejection_status_vector[2]+1,rejection_status_vector[3]+1,null_hyp,stage]){
				stage_when_hypotheses_rejected[null_hyp] <- stage;
				rejection_status_vector[null_hyp] <- 1;
				new_null_hypothesis_rejected_current_stage <- TRUE;
				break;
			}
		     }
		}
	    }	   				 
	    ## Enrollment Modification Rule
	    if(enrollment_status[1]==1 && (stage_when_hypotheses_rejected[1]==stage || Z_subpop_1_cumulative[stage,simulated_trial_iteration] < H01_futility_boundaries[stage] || Z_combined_population_cumulative[stage,simulated_trial_iteration] < H0C_futility_boundaries[stage])){ # stop enrollment from subpop 1
			final_stage_subpop_1_enrolled_up_through[simulated_trial_iteration] <- stage;
			enrollment_status[1] <- 0;					   
	    }
	    if(enrollment_status[2]==1 && (stage_when_hypotheses_rejected[2]==stage || Z_subpop_2_cumulative[stage,simulated_trial_iteration] < H02_futility_boundaries[stage] || Z_combined_population_cumulative[stage,simulated_trial_iteration] < H0C_futility_boundaries[stage])){ # stop enrollment from subpop 2
			final_stage_subpop_2_enrolled_up_through[simulated_trial_iteration] <- stage;
			enrollment_status[2] <- 0;					   
	    }
	    if(enrollment_status[3]==1 && (enrollment_status[1]==0 || enrollment_status[2]==0)){
			enrollment_status[3] <- 0;					   
	    }		
	}
	reject_H01[simulated_trial_iteration] <- rejection_status_vector[1];
	reject_H02[simulated_trial_iteration] <- rejection_status_vector[2];
	reject_H0C[simulated_trial_iteration] <- rejection_status_vector[3];
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
'H01_efficacy_boundaries'=all_efficacy_boundaries[1,1,1,1,], # efficacy boundaries for null hypothesis H01
'H02_efficacy_boundaries'=all_efficacy_boundaries[1,1,1,2,], # efficacy boundaries for null hypothesis H02
'H0C_efficacy_boundaries'=all_efficacy_boundaries[1,1,1,3,], # efficacy boundaries for null hypothesis H0C
'all_efficacy_boundaries'=all_efficacy_boundaries,
'H01_futility_boundaries'=H01_futility_boundaries, # futility boundaries for null hypothesis H01
'H02_futility_boundaries'=H02_futility_boundaries, # futility boundaries for null hypothesis H02
'H0C_futility_boundaries'=H0C_futility_boundaries, # futility boundaries for null hypothesis H0C
SS_dist = SS_dist, #Full distribution of SS
dur_dist = dur_dist #Full distribution of trial duration
))
	 
}




