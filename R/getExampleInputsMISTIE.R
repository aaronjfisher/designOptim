

#' Example inputs for \code{\link{optimizeTrial}}
#'
#' @param max_K maximum number of stages
#' @param trial_method either 'cov' or 'MB'
#' @param bound_style either 'unstructured' or 'structured'
#' @param which_obj either 'ss' or 'dur' to minimize expected sample size or duration respectively
#' @param iter number of search iterations 
#' @param search_stage1 whether to pre-calculate an optimal 1-stage trial
#' @param ... passed to \code{\link{min_n_feasible}}
#' @export
#' @seealso \code{\link{optimizeTrial}}
getExampleInputsMISTIE <-function(
	max_K = 4, 
	trial_method='cov', 
	bound_style='unstructured', #
	which_obj = 'ss',
	iter=50,
	search_stage1 = FALSE,
	...
	){


	num_stages = max_K

	##### MISTIE parameters

	p1 <- 1/3 #subpopulation 1 proportion
	delta <- 0.122 #treatment effect

	mean_s1_con <- 0.29
	mean_s2_con <- 0.29
	mean_s1_trt <- mean_s1_con + delta
	mean_s2_trt <- mean_s2_con + delta

	var_s1_con <- mean_s1_con * (1-mean_s1_con)
	var_s2_con <- mean_s2_con * (1-mean_s2_con)
	var_s1_trt <- mean_s1_trt * (1-mean_s1_trt)
	var_s2_trt <- mean_s2_trt * (1-mean_s2_trt)

	enroll_MISTIE = 420
	delay_MISTIE = 180/365 # (in years)

	##### 



	################################
	# Define prior over which to calculate the power & expected sample size.
	# We set the prior on (\delta_1, \delta_2) to be at
	# point masses. These points are labelled: a1a2, a1n2, 
	# n1a2, n1n2. "aX" denotes that the alternative 
	# hypothesis holds for subgroup "X," and "nX" denotes 
	# that the null hypothesis holds for subgroup "X".


	cases <- list()

	power_min_const <- 0.8 #Desired power

	# We consider four cases
	cases$a1a2 <-
	cases$a1n2 <-
	cases$n1a2 <- 
	cases$n1n2 <- list()

	for( i in 1:length(cases)){

		mins_list <-list()
		mins_list$Pow_H0C <- 0 #Initial values that will soon be overriden
		mins_list$Pow_H01 <- 0
		mins_list$Pow_H02 <- 0

		cases[[i]]$power_mins <- mins_list


		cases[[i]]$var_s1_trt <- var_s1_trt
		cases[[i]]$var_s1_con <- var_s1_con
		cases[[i]]$var_s2_trt <- var_s2_trt
		cases[[i]]$var_s2_con <- var_s2_con

		cases[[i]]$mean_s1_con <- mean_s1_con
		cases[[i]]$mean_s2_con <- mean_s2_con

		cases[[i]]$weight <- 1/length(cases)
	}

	#required power for each hypothesis in each case
	cases$a1a2$power_mins$Pow_H0C <- power_min_const 
	cases$a1n2$power_mins$Pow_H01 <- power_min_const
	cases$n1a2$power_mins$Pow_H02 <- power_min_const

	cases$a1a2$mean_s1_trt <-
	cases$a1n2$mean_s1_trt <- mean_s1_trt

	cases$a1a2$mean_s2_trt <-
	cases$n1a2$mean_s2_trt <- mean_s2_trt

	cases$n1n2$mean_s2_trt <- 
	cases$a1n2$mean_s2_trt <- mean_s2_con
	cases$n1n2$mean_s1_trt <-
	cases$n1a2$mean_s1_trt <- mean_s1_con


	# str(cases)


	################################
	# Search for the minimum feasible sample size


	#For "MB" method we need to specify graph parameters
	stage1list <-list()
	if(trial_method!='cov'){
		stage1list <-list( graph_edge_12=1/2,graph_edge_2C=1/2, graph_edge_C1=1/2)
	}

	if(!search_stage1){
		stage1_feasible <- NULL
		n_total <- 2100
	}else{
		cat('Searching for smallest 1-stage trial that meets constraints...\n')
		(stage1_feasible <-min_n_feasible(
				FWER=0.025,
				p1=p1,
				cases=cases,
				min_n=500,
				max_n=2000,
				step_n=20, 
				trial_method=trial_method,
				trial_args=stage1list,
				...
				# npoints_sqrt=10, #Set artifically low for example purposes
				# showiter=showiter,
			))


		n_feasible_multiplier <- 1.25
		n_total <- stage1_feasible['n_total']  * n_feasible_multiplier
		names(n_total)<-c()
	}

	n_per_stage <- rep(n_total/num_stages, max_K)
		#only the first `num_stages` elements of this vector will be used, but it must be of length = max_K.



	##############################
	# Set up argument list to send to optimizeTrial

	# `args_list` contains arguments to be passed to buildTrial. Some
	# of these are fixed, and some are to be optimized over. For those
	# variables to optimize over, we enter the initial/starting 
	# values into `args_list`.

	# IMPORTANT NOTE: All vectors must be of length max_K, 
		# even if the initial num_stages is < max_K.
		# only the first K' elements will be used, where K' is the proposed 	
		# value of num_stages at any iteration of the SANN search.
		# num_stages must be a positive integer <= max_K.
		# Setting num_stages=1 will yeild a non-adaptive design.


	args_list<-list(
		#All elements must be numeric vectors (or scalars)
		#All per-stage vectors *must* have length = max_K (see note above)

		num_stages = num_stages, 	
	  	n_per_stage = n_per_stage,
		n_total = n_total,

		p1 = p1, 
		r1 = 1/2,
		r2 = 1/2,

		FWER = 0.025,
		enrollment_rate_combined = enroll_MISTIE,
		delay = delay_MISTIE,

		graph_edge_12 = 1/2, #alpha reallocation graph -- not used in "cov" design.
		graph_edge_2C = 1/2,
		graph_edge_C1 = 1/2,

		time_limit = 90, #seconds

		errtol = .01,
		#following two arguments for GanzBretz are overridden by build_precision option
		abseps = 0.000001, 
		maxpts = 10000,

		iter = iter  

		#Note: FWER_allocation_matrix should not be included, as this is not a vector.
	)

	if(bound_style=='structured'){
		#here, the parameter length of 3 is hardcoded to match H01, H02 & H0C (not to match max_K)
		args_list$delta_futility <- c(0.5,0.5,0.5)
		args_list$intercepts_futility <- c(0,0,0)
		args_list$H01_futility_boundary_const <- -1
		args_list$H02_futility_boundary_const <- -1
		args_list$H0C_futility_boundary_const <- -2

		args_list$delta_eff <-  c(0.5,0.5,0.5)
		args_list$H01_eff_total_allocated <- 0.5 #Note, setting any of these = 1 here will cause errors if we also use logit_search here
		args_list$H02_eff_total_allocated <- 0.5
		args_list$H0C_eff_total_allocated <- 0.5
	}
	if(bound_style=='unstructured'){
		args_list$H01_eff_allocated <- rep(.5,max_K)
		args_list$H02_eff_allocated <- rep(.5,max_K)
		args_list$H0C_eff_allocated <- rep(.5,max_K)

		args_list$H01_futility_boundaries <- rep(0,max_K)
		args_list$H02_futility_boundaries <- rep(0,max_K)
		args_list$H0C_futility_boundaries <- rep(-4,max_K)
	}



	# Vector that tells which arguments to optimize over:
	index2optim <- c(
		'n_total',
		'n_per_stage',
		'H01_eff_allocated',
		'H02_eff_allocated',
		'H0C_eff_allocated',
		'graph_edge_12',
		'graph_edge_2C',
		'graph_edge_C1',
		'H01_futility_boundaries',
		'H02_futility_boundaries',
		'H0C_futility_boundaries',
		'delta_futility',
		'intercepts_futility',
		'H01_futility_boundary_const',
		'H02_futility_boundary_const',
		'H0C_futility_boundary_const',
		'H01_eff_total_allocated',
		'H02_eff_total_allocated',
		'H0C_eff_total_allocated',
		'delta_eff'
	) # Note: Quantities describing simulation precision
	  # (e.g. abseps, maxpts) should never be put here.

	index2optim <- intersect(names(args_list),index2optim)
	index2fix <- names(args_list)[!names(args_list)%in%index2optim]

	args_list_init <- args_list[index2optim]
	args_list_fixed <- args_list[index2fix]


	# Arguments to search for on the logit space
	logit_search<-c(
		"H01_eff_allocated",
		"H02_eff_allocated",
		"H0C_eff_allocated",
		"H01_eff_total_allocated",
		"H02_eff_total_allocated",
		"H0C_eff_total_allocated",
		'graph_edge_12',
		'graph_edge_2C',
		'graph_edge_C1'
	)


	if(which_obj == 'ss') obj_fun <- min_E_SS_power_constraints
	if(which_obj ==  'dur') obj_fun <- min_E_dur_power_constraints

	return(list(
		'args_list_init' = args_list_init,
		'args_list_fixed' = args_list_fixed,
		'cases' = cases,
		'local_n_search' = FALSE,
		'objective_fun' = obj_fun,  
		'max_K' = max_K,
		'trial_method'= trial_method,
		'optim_method' = 'SANN',
		'maxit' = iter,
		'logit_search' = logit_search,
		'parscale_ratio_N' = 100,
		'parscale_ratio_K' = 1,
		'parscale_ratio_logit' = 3,
		'print_interval' = ceiling(iter/10),
		'returnFullPath' = TRUE,
		'stage1_feasible'=stage1_feasible,
		'verbose' = TRUE,
		'print_obj' = TRUE,
		'build_precision'=FALSE #maintain high precision throughout search
		))
}



