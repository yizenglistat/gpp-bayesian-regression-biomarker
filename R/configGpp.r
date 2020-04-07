gpp_config_fun<-function(t, nknots=100, nbeta, phi_sd, kappa=2, step=0.01){
	N 				<- length(t)						# number of individuals
	eta_idx 		<- rep(NA, N)						# linear predictor eta initialized

	t_lower 		<- min(t)							# lower bound of variable t
	t_upper 		<- max(t)							# upper bound of variable t
	
	t_unique 		<- unique(sort(t))					# sorted extracted t_unique
	nunique 		<- length(t_unique)					# number of unique in t sequence 
	
	t_new   		<- seq(t_lower, t_upper, by=step)	# new t sequence
	nnew 			<- length(t_new)					# number of unique for new t sequence
														# generate indicator matrix for tunqiue for N individuals 
	for(idx in 1:nunique){
		eta_idx[ t == t_unique[idx] ] = idx 			# store those individuals for each t_unique
	}
	E_mat = Matrix(0, N, nunique)						# initialize E_mat indicator matrix N x nunique
	E_mat[cbind(1:N, eta_idx)] = 1 						# assign those individuals to 1 if appears in eta_idx

	# distance matrix for computing correlation matrix later in convenience
	t_knots = seq(t_lower, t_upper, length.out=nknots)
	unique_dist = as.matrix(dist(t_unique, t_unique, method = "manhattan", diag = FALSE, upper = FALSE))
	knots_dist  = as.matrix(dist(t_knots, t_knots, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist  = matrix(NA, nunique, nknots)
	cross_dist_new = matrix(NA, nnew, nknots)
	for(idx in 1:nknots){
	  cross_dist[,idx] = abs(t_unique-t_knots[idx])
	  cross_dist_new[,idx] = abs(t_new-t_knots[idx])
	}

	# configure the range of phi
	t_range = t_upper-t_lower
	phi_seq = seq(0.01,100,by=0.01)
	phi_upper = phi_lower = rep(NA,length(phi_seq))
	for(idx in 1:length(phi_seq)){
	  phi_lower[idx] = matern(t_range*2/30,phi=phi_seq[idx],kappa = kappa)
	  phi_upper[idx] = matern(t_range*2/3,phi=phi_seq[idx],kappa = kappa)
	}
	phi_min = phi_seq[which.min(abs(phi_lower-0.05))]
	phi_max = phi_seq[which.min(abs(phi_upper-0.05))]

	# configure required beta functional parameters 

			# ---------------------- beta_config matrix description ----------------------- #
			# each columns means a beta function											#
			# row1: gamma distribution shape, prior distribution for tau 					#
			# row2: gamma distribution scale, prior distribution for tau 					#
			# row3: sampled tau from gamma prior distribution based on row1 and row2 		#
			# row4: uniform distribution lower bound for phi 								#
			# row5: uniform distribution upper bound for phi 								#
			# row6: phi averaged value based on uniform distribution based on row4 and row5	#
			# row7: transformed phi for M-H sampler, which made a transformation on row6 	#
			# row8: hyperparameter phi_sd, default 0.1										#
			# row9: hyperparameter kappa, default 2											#
			# ----------------------------------------------------------------------------- #

	beta_config <- matrix(NA, 9, nbeta)
	for (d in 1:nbeta){
		beta_config[1,d] <- 2
		beta_config[2,d] <- 1
		beta_config[3,d] <- rgamma(1, beta_config[1,d], beta_config[1,d])
		beta_config[4,d] <- phi_min
		beta_config[5,d] <- phi_max
		beta_config[6,d] <- (beta_config[4,d]+beta_config[5,d])/2
		beta_config[7,d] <- logit((beta_config[6,d]-beta_config[4,d])/(beta_config[5,d]-beta_config[4,d]))
		beta_config[8,d] <- phi_sd
		beta_config[9,d] <- kappa
	}

	R_knots_list 	 	<- list()						# initial list knots corr mat
	R_knots_inv_list	<- list() 						# initial list inverse knots corr mat
	R_cross_list 		<- list()						# initial list nunique x nknots cross corr mat
	Q_list 				<- list()						# initial list conditional transform mat
	for (d in 1:nbeta){
		R_knots_list[[d]]		<- toeplitz(matern(knots_dist[1,], phi= beta_config[6,d], kappa = kappa))
		R_knots_inv_list[[d]] 	<- TrenchInverse(R_knots_list[[d]])
		R_cross_list[[d]] 		<- matern(cross_dist, phi=beta_config[6,d], kappa=kappa)
		Q_list[[d]] 			<- R_cross_list[[d]]%*%R_knots_inv_list[[d]]
	}

	# output list
	output_list = list(E_mat=E_mat,kappa=kappa, t_new=t_new, nnew=nnew,
					   t_unique=t_unique, nunique=nunique, t_knots=t_knots, nknots=nknots,
					   unique_dist=unique_dist, 
					   knots_dist=knots_dist, 
					   cross_dist=cross_dist,
					   cross_dist_new=cross_dist_new, 
					   phi_init=c(phi_min,phi_max),
					   beta_config=beta_config,
					   R_knots_list=R_knots_list, 
					   R_knots_inv_list=R_knots_inv_list,
					   R_cross_list=R_cross_list,
					   Q_list=Q_list
					   )
	return(output_list)
}