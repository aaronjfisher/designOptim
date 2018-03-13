

#
#' A simulated annealing variant that re-evaluates the reference point at each comparison
#'
#' @param fn function to minimize
#' @param par initial parameters
#' @param  control analogous to \code{\link{optim}}
#' @export
#' @import mvtnorm
#' @importFrom stats rbinom
#' @return optimized value, found at minimum of all evaluated solutions.
#' @examples \dontrun{
#' fn<-function(x){
#'     abs(x) + rnorm(n=1,mean=0,sd=.3) 
#' }
#' par = 2
#' control <- list()
#' control$parscale <- 2
#' control$maxit <- 3000
#' control$temp<-10
#' sann_random(fn, par, control)
#'}
sann_random <- function(fn,par, control=list()){
	
	if(is.null(control$parscale)) control$parscale<-rep(1,length(par))
	if(is.null(control$maxit)) control$maxit<-10000
	if(is.null(control$temp)) control$temp<-10000
	temp_start <- control$temp

	theta_c <- par #current theta
	theta_log <- list()
	fun_log <- rep(NA,control$maxit)

	for(i in 1:control$maxit){
		#Temperature t
		t <- temp_start / log((i-1)  + exp(1)) #this is equivalent to setting tmax=1 in optim 

		#proposal
		sigma <- control$parscale*(t/temp_start)
		if(length(control$parscale)>1) sigma <- diag(control$parscale)*(t/temp_start)
		theta_p <- rmvnorm(1, mean=theta_c, sigma=as.matrix(sigma)) #proposal
			#this sigma is not fully explained in ?optim,
			#(a proportionality constant is not given)
			#This form of sigma is motivated by trying
			#to mimic the proposal pattern of optim.

		#Core aspect of sann_random: we evaluate objective function twice
		fp<-fn(theta_p)
		fc<-fn(theta_c)

		theta_log <- c(theta_log, theta_p)
		fun_log[i] <- fp

		prob_accept <- exp( - (fp - fc) / t )

		acpt <- rbinom(1, 1, prob=min(1,prob_accept))==1
		if(acpt) theta_c <- theta_p

	}

	# plot(unlist(theta_log),type='l')
	# plot(fun_log,type='l')


	theta_log[[which(fun_log==min(fun_log))]]
		#this is especially prone to the monte-carlo error problem.
		#but works well here better than returning final value of theta_c
}