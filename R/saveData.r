saveData<-function(MC_list, taskid=character(0), reps=reps, rep=1, visual=FALSE, savepath='./output/'){
	# -$-$-$-$-$-$-$-$- posterior inference on dorfman test data -$-$-$-$-$-$-$-$-$-$-$-$-$-$-$

	t_knots 					<- MC_list$t_knots
	beta_MC 					<- MC_list$beta 							# Monte Carlo samples, list with length nbeta and each element of list is N x nkeep matrix
	beta_new_MC 				<- MC_list$beta_new 						# Monte Carlo samples, list with length nbeta and each element of list is nnew x nkeep matrix
	beta_knots_MC 				<- MC_list$beta_knots 						# Monte Carlo samples, list with length nbeta and each element of list is nnew x nkeep matrix
	tau_MC 						<- MC_list$tau 								# Monte Carlo samples, list with length nbeta and each element of list is 1 x nkeep matrix
	phi_MC 						<- MC_list$phi 								# Monte Carlo samples, Monte Carlo samples, list with length nbeta and each element of list is 1 x nkeep matrix
	Y_MC						<- MC_list$Y
	sigma_MC 					<- MC_list$sigma

	N 							<- nrow(beta_MC[[1]])
	nbeta 						<- length(beta_MC)							# number of beta functions
	nnew 						<- length(MC_list$t_new)					# length of t_new
	nknots 						<- nrow(MC_list$beta_knots[[1]])


	theta 						<- sigma_MC

	
	for(d in 1:nbeta){
		theta 					<- rbind(theta, tau_MC[[d]])
	}

	for(d in 1:nbeta){
		theta 					<- rbind(theta, phi_MC[[d]])
	}

	beta 						<- c()
	beta_new 					<- c()
	beta_knots 					<- c()
	for(d in 1:nbeta){
		beta 					<- rbind(beta, beta_MC[[d]])
		beta_new 				<- rbind(beta_new, beta_new_MC[[d]])
		beta_knots				<- rbind(beta_knots, beta_knots_MC[[d]])
	}

	# # -$-$-$- posterior mean on dorfman test data -$-$-$
	theta_mean 					<- apply(theta, 1, mean)							# posterior mean for Se, Sp, tau and phi, a vector of length (2*nassay+2)
	beta_mean 					<- matrix(apply(beta,1,mean),N,nbeta)				# posterior mean for beta functional, a matrix of N x nbeta
	beta_new_mean 				<- matrix(apply(beta_new,1,mean),nnew,nbeta)		# posterior mean for beta_new functional, a matrix of nnew x nbeta
	beta_knots_mean 			<- matrix(apply(beta_knots,1,mean),nknots,nbeta)
	# # -$-$-$- posterior median on dorfman test data -$-$-$
	theta_med 					<- apply(theta, 1, median)							# posterior median for Se, Sp, tau and phi, a vector of length (2*nassay+2)
	beta_med 					<- matrix(apply(beta,1,median),N,nbeta)				# posterior median for beta functional, a matrix of N x nbeta
	beta_new_med				<- matrix(apply(beta_new,1,median),nnew,nbeta)		# posterior median for beta_new functional, a matrix of nnew x nbeta
	beta_knots_med				<- matrix(apply(beta_knots,1,median),nknots,nbeta)

	# # -$-$-$- posterior standard error on dorfman test data -$-$-$
	theta_sd 					<- apply(theta, 1, sd)								# posterior standard error for Se, Sp, tau and phi, a vector of length (2*nassay+2)
	beta_sd 					<- matrix(apply(beta,1,sd),N,nbeta)					# posterior standard error for beta(t) functional, a matrix of N x nbeta
	beta_new_sd  				<- matrix(apply(beta_new,1,sd),nnew,nbeta)			# posterior standard error for beta(t_new) functional, a matrix of nnew x nbeta
	beta_knots_sd  				<- matrix(apply(beta_knots,1,sd),nknots,nbeta)
	# # -$-$-$- posterior credible interval on dorfman test data -$-$-$
	theta_ci 					<- HPDinterval(as.mcmc(t(theta)))					# posterior credible interval for Se, Sp, tau and phi, a vector of length (2*nassay+2)
	beta_ci 					<- HPDinterval(as.mcmc(t(beta)))					# posterior credible interval for beta functional, a matrix of N x nbeta
	beta_new_ci 				<- HPDinterval(as.mcmc(t(beta_new)))				# posterior credible interval for beta_new functional, a matrix of nnew x nbeta
	beta_knots_ci 				<- HPDinterval(as.mcmc(t(beta_knots)))
	# # summary all posterior inference
	beta_ci_lower 				<- matrix(beta_ci[,1],N,nbeta)						# posetrior credible interval lower bound for beta functional, a matrix of N x nbeta
	beta_ci_upper				<- matrix(beta_ci[,2],N,nbeta)						# posetrior credible interval upper bound for beta functional, a matrix of N x nbeta
	beta_new_ci_lower 			<- matrix(beta_new_ci[,1],nnew,nbeta)				# posetrior credible interval lower bound for beta_new functional, a matrix of nnew x nbeta
	beta_new_ci_upper			<- matrix(beta_new_ci[,2],nnew,nbeta)				# posetrior credible interval upper bound for beta_new functional, a matrix of nnew x nbeta
	beta_knots_ci_lower 		<- matrix(beta_knots_ci[,1],nknots,nbeta)				# posetrior credible interval lower bound for beta_new functional, a matrix of nnew x nbeta
	beta_knots_ci_upper			<- matrix(beta_knots_ci[,2],nknots,nbeta)

	theta_summ 					<- data.frame(cbind(theta_mean, theta_med, theta_sd, theta_ci))
	colnames(theta_summ)		<- c('mean', 'median', 'sd', 'lower', 'upper')

	beta_summ 					<- data.frame(cbind(beta_mean, beta_med, beta_sd, beta_ci_lower, beta_ci_upper))
	colnames(beta_summ)			<- c(paste0('mean_beta',0:(nbeta-1)), paste0('med_beta',0:(nbeta-1)),
									paste0('sd_beta',0:(nbeta-1)),paste0('lower_beta',0:(nbeta-1)),paste0('upper_beta',0:(nbeta-1)))

	beta_new_summ 				<- data.frame(cbind(beta_new_mean, beta_new_med, beta_new_sd, beta_new_ci_lower, beta_new_ci_upper))
	colnames(beta_new_summ)		<- c(paste0('mean_beta',0:(nbeta-1)), paste0('med_beta',0:(nbeta-1)),
									paste0('sd_beta',0:(nbeta-1)),paste0('lower_beta',0:(nbeta-1)),paste0('upper_beta',0:(nbeta-1)))

	beta_knots_summ 				<- data.frame(cbind(t_knots,beta_knots_mean, beta_knots_med, beta_knots_sd, beta_knots_ci_lower, beta_knots_ci_upper))
	colnames(beta_knots_summ)		<- c('t_knots',paste0('mean_beta',0:(nbeta-1)), paste0('med_beta',0:(nbeta-1)),
									paste0('sd_beta',0:(nbeta-1)),paste0('lower_beta',0:(nbeta-1)),paste0('upper_beta',0:(nbeta-1)))

	postfix <- reps*(taskid-1)+rep
	fwrite(theta_summ, file = paste0(savepath,'theta_summ',paste0(rep('0',3-nchar(postfix)),collapse = ''),postfix,'.csv'), row.names=FALSE)
	#write.csv(beta_summ, file = paste0(savepath,'beta_summ',reps*(taskid-1)+rep,'.csv'), row.names=FALSE)
	#write.csv(beta_new_summ, file = paste0(savepath, 'beta_new_summ',reps*(taskid-1)+rep,'.csv'), row.names=FALSE)
	fwrite(beta_knots_summ, file = paste0(savepath,'beta_knots_summ',paste0(rep('0',3-nchar(postfix)),collapse = ''),postfix,'.csv'), row.names=FALSE)

}

