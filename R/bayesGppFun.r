bayesGppFun <- function(center=FALSE, Xd, omega=1, sigma=1, h, h_beta, E_mat, Q_mat, knots_dist, cross_dist, R_knots_mat, R_knots_inv_mat, R_cross_mat, beta_d_config){
	
	# ------------------------------------------ Update beta_d(t), beta_d(t_unique), beta_d(t_knots) ------------------------------------------ #

	tau  			 <- beta_d_config[3]													# tau parameter for beta_d
	Omega_d2		 <- Diagonal(x=Xd^2*omega)												# Omega_d2. diagonal matrix
	EOmegah_beta  	 <- crossprod(E_mat, Xd*omega*h_beta)									# part of updated mean of GPP
	D_mat 			 <- crossprod(E_mat,Omega_d2)%*%E_mat
	presicion_mat 	 <- forceSymmetric(crossprod(Q_mat, D_mat)%*%Q_mat + R_knots_inv_mat*tau*(sigma**2))
	presicion_decomp <- chol(presicion_mat)													# chol decomposition to make inverse faster
	Sigma_beta_d  	 <- chol2inv(presicion_decomp) 											# updated covariance matrix of GPP
	mu_beta_d  		 <- Sigma_beta_d %*% crossprod(Q_mat,EOmegah_beta)   					# updated mean of GPP

	# use override old matrix to make name short, otherwise, updated_beta_d_knots, etc
	nknots 		  	 <- nrow(R_knots_mat)													# number of knots
	beta_d_knots  	 <- as.numeric(backsolve(presicion_decomp/sigma,rnorm(nknots))+mu_beta_d)		# sample from presicion matrix,solve upper triangle system
	beta_d_unique 	 <- Q_mat%*%beta_d_knots 												# convert from knots to unique values in t
	beta_d_knots  	 <- beta_d_knots  - mean(beta_d_unique)*center							# center beta_d(t_knots) or not 
	beta_d_unique 	 <- beta_d_unique  - mean(beta_d_unique)*center							# center beta_d(t_unique) or not
	beta_d        	 <- as.vector(E_mat%*%beta_d_unique)									# updated beta_d(t)

	# -------------------------------------------------------------- Update tau -------------------------------------------------------------- #
	
	a 				 <- beta_d_config[1]													# shape parameter for gamma prior
	b 				 <- beta_d_config[2]													# scale parameter for gamma prior
	updated_a 		 <- a + nknots/2														# updated shape and scale for gamma prior
	updated_b 		 <- as.numeric(b + crossprod(beta_d_knots, R_knots_inv_mat)%*%beta_d_knots/2)
	updated_tau 	 <- rgamma(1, updated_a, updated_b)										# sample a new tau from updated gamma prior
	
	# -------------------------------------------------------- Update GPP correlation matrix ------------------------------------------------------- #
	
	# -$-$-$-$-$-$-$-$- configure required arguments -$-$-$-$-$-$-$-$-

	phi_min 							<- beta_d_config[4]									# phi lower bound
	phi_max 		 					<- beta_d_config[5]									# phi upper bound
	phi     						 	<- beta_d_config[6]									# phi parameter in covariance matrix
	trphi   						  	<- beta_d_config[7]									# trphi value for M-H algorithm.
	phi_sd 								<- beta_d_config[8]									# hyperparameter, phi_sd
	kappa 								<- beta_d_config[9]									# hyperparameter, kappa

	# -$-$-$-$-$-$-$-$- proposed point from simple random walk -$-$-$-$-$-$-$-$-

	proposed_trphi 						<- rnorm(1, trphi, phi_sd)										# propose a new phi from a Gaussian random walk		
	proposed_phi   						<- (exp(proposed_trphi)*phi_max+phi_min)/(1+exp(proposed_trphi)) 	# transformation on the proposed phi
	proposed_R_knots_mat 				<- toeplitz(matern(knots_dist[1,], phi=proposed_phi, kappa=kappa)) 	# updated knots mat for beta_d
	proposed_R_knots_inv_mat 			<- TrenchInverse(proposed_R_knots_mat)								# updated inverse knots mat for beta_d
	proposed_R_cross_mat 				<- matern(cross_dist, phi=proposed_phi, kappa=kappa)				# updated nunique x nknots cross mat for beta_d
	proposed_Q_mat 						<- proposed_R_cross_mat%*%proposed_R_knots_inv_mat					# updated Q_mat for beta_d

	h_minus_eta 						<- h_beta - Xd * beta_d 												# h - eta where eta is linear predictor
	proposed_h_minus_eta 				<- h_beta - Xd * (E_mat%*%proposed_Q_mat%*%beta_d_knots)				# proposed (h - eta): N x 1 vector

	proposed_one 						<-  -determinant(proposed_R_knots_mat)$modulus/2
		  									-updated_tau*crossprod(beta_d_knots, proposed_R_knots_inv_mat)%*%beta_d_knots/2
		  									-crossprod(proposed_h_minus_eta*omega, proposed_h_minus_eta)/(2*sigma**2)
	
	current_one 						<-  -determinant(R_knots_mat)$modulus/2
									  		-updated_tau*crossprod(beta_d_knots, R_knots_inv_mat)%*%beta_d_knots/2
									  		-crossprod(h_minus_eta*omega, h_minus_eta)/(2*sigma**2)

	rate 								<- (proposed_phi-phi_min)*(phi_max-proposed_phi)/((phi_max-phi)*(phi-phi_min))	# to control the accept probability								
	accept_prob 						<- min(1, as.numeric(exp(proposed_one-current_one)*rate))						# accept probability to proposed one or reject

	bernoulli_var 						<- rbinom(1,1,accept_prob)														# generate bernoulli variable with accept prob
	updated_phi  						<- bernoulli_var*proposed_phi   + (1-bernoulli_var)*phi 						# updated phi (could be no change)
	updated_trphi						<- bernoulli_var*proposed_trphi + (1-bernoulli_var)*trphi 						# updated trphi (could be no change)

	# use override old matrix to be consistent with the input matrix name
	if(bernoulli_var==1){ 																	# if proposed accepted, then update all matrics
	  R_knots_mat 						<- proposed_R_knots_mat 							# updated knots correlation matrix
	  R_knots_inv_mat 					<- proposed_R_knots_inv_mat							# updated inverse knots correlation matrix			
	  R_cross_mat 						<- proposed_R_cross_mat 							# updated (nunique x nknots) cross unique knots correlation matrix
	  Q_mat 							<- proposed_Q_mat 									# updated Q_mat converting knots to unique values in t.
	}

	output_list <- list(beta_d=beta_d, beta_d_knots=beta_d_knots, beta_d_unique=beta_d_unique, 
						tau=updated_tau, phi=updated_phi, trphi=updated_trphi,
						R_knots_mat=R_knots_mat, R_knots_inv_mat=R_knots_inv_mat, R_cross_mat=R_cross_mat,Q_mat=Q_mat)

	return(output_list)
}