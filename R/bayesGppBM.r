bayesGppBM <- function(X, G, t, Y_pool, Z, I, W, size, visual=FALSE,
					true_params=list(Y_indv=Y_indv,sigma_true=sigma_true,beta_true=beta_true,alpha_true=alpha_true), # this is just for plotting
					hyper_params=list(phi_sd=0.1, kappa=2,nknots=50,sigma_init=1),
					gibbs_config=list(nburn=1000,nkeep=2000,nthin=5,nstep=50)){

	# ============================================================================================= #
	# ...................................... CONFIGURE BLOCK ...................................... #
	# ============================================================================================= #
	# cat("\014")										# clear the screen
	# --------------------------------- configure data information -------------------------------- #
	
	N	  	<- nrow(X) 									# number of individuals/patients
	npools 	<- nrow(Z) 									# number of pools
	cj 		<- rep(size,npools)
	nbeta 	<- ncol(X)									# number of unknown beta fun including beta0
	nalpha 	<- ncol(G) 									# number of linear effect coefficients

	# --------------------------------- configure true parameters -------------------------------- #

	Y_indv 			<- true_params$Y_indv
	sigma_true 		<- true_params$sigma_true
	beta_true 		<- true_params$beta_true
	alpha_true 		<- true_params$alpha_true

	# --------------------------------- configure Gibbs sampler --------------------------------- #

	nburn 	<- gibbs_config$nburn						# number of burn in Gibbs sampling
	nkeep 	<- gibbs_config$nkeep						# number of keep in Gibbs sampling
	nthin 	<- gibbs_config$nthin						# number of thinning in Gibbs sampling
	nstep 	<- gibbs_config$nstep						# control visualization in Gibbs sampling

	# --------------------------------------- configure GPP -------------------------------------- #

	phi_sd 			<- hyper_params[[1]]				# phi_sd control smooth about corr mat
	kappa			<- hyper_params[[2]]				# kappa control smooth about corr mat
	nknots			<- hyper_params[[3]]

	gpp_config 		<- gpp_config_fun(t=t, 				# t sequence
								 	  nknots=nknots, 	# number of selected knots
								 	  nbeta=nbeta,		# number of beta including intercept
								 	  phi_sd=phi_sd,	# hyperparameter phi_sd, default is 0.1
								 	  kappa=kappa)		# hyperparameter kappa, default is 2

	t_new 			<- gpp_config$t_new					# t_new seq, all unique values inside domain t
	nnew 			<- gpp_config$nnew 					# *note* number of t_new = number of t_unique 

	t_unique 		<- gpp_config$t_unique 				# unique t values in t sequence
	nunique 		<- gpp_config$nunique 				# number of unique t values
	t_knots  		<- gpp_config$t_knots 				# initial selected knots
	nknots 			<- gpp_config$nknots				# number of selected knots
	
	# in order to calculate correlation matrix of GPP easily, compute distance matrix in advance.
	unique_dist 	<- gpp_config$unique_dist			# nunique x nunique, L1 dist matrix
	knots_dist  	<- gpp_config$knots_dist			# nknot x nknot, L1 dist matrix
	cross_dist  	<- gpp_config$cross_dist 			# nunique x nknot, L1 dist matrix
	cross_dist_new 	<- gpp_config$cross_dist_new 		# nnew x nknot, L1 dist matrix

	phi_init 		<- gpp_config$phi_init				# uniform prior, uniform(phi_min, phi_max)
	phi_min  		<- phi_init[1]						# uniform lower bound
	phi_max  		<- phi_init[2]						# uniform upper bound
	
			# ---------------------- beta_config matrix description ----------------------- #
			# each columns means a beta function											#
			# row1: gamma distribution shape, prior distribution for tau 					#
			# row2: gamma distribution scale, prior distribution for tau 					#
			# row3: sampled tau from gamma prior distribution based on row1 and row2 		#
			# row4: uniform distribution lower bound for phi 								#
			# row5: uniform distribution upper bound for phi 								#
			# row6: phi averaged value based on uniform distribution based on row4 and row5	#
			# row7: proposed phi for M-H sampler, which made a transformation on row6 		#
			# row8: hyperparameter phi_sd, default 0.1										#
			# row9: hyperparameter kappa, default 2											#
			# ----------------------------------------------------------------------------- #

	beta_config 	<- gpp_config$beta_config 			# matrix of beta configuration for gpp 

	R_knots_list 	<- gpp_config$R_knots_list 			# list of correlation matrix, nknots x nknots
	R_knots_inv_list<- gpp_config$R_knots_inv_list		# list of inverse correlation matrix, nknots x nknots
	R_cross_list 	<- gpp_config$R_cross_list 			# list of cross correlation matrix, nunique x nknots
	Q_list 			<- gpp_config$Q_list 				# convert beta(t_knots) vector back to beta(t_unique) vector

	beta			<- matrix(0, N, nbeta)				# initial beta(t) matrix
	beta_unique 	<- matrix(0, nunique, nbeta)		# initial beta(t_unique) matrix
	beta_knots 	 	<- matrix(0, nknots, nbeta)			# initial beta(t_knots) matrix
	Xbeta 			<- rep(0,N)
	Xbeta 			<- createXbeta(N,nbeta,Xbeta,X,beta)

	# ----------------------------------- configure shift alpha ------------------------------------- #
	
	alpha 			<- c(1,rep(0,nalpha-1))
	Galpha 			<- G%*%alpha

	# ----------------------------------- configure Latent Y ------------------------------------- #
	  
	sigma   		<- hyper_params[[4]]
	Y    			<- rep(Y_pool, cj)

	# ============================================================================================= #
	# .................................... GIBBS SAMPLING BLOCK ................................... #
	# ============================================================================================= #

	#------------------------------------- Start Gibbs sampling ------------------------------------#
	
	kept_beta			<- lapply(1:nbeta, matrix, data = NA, nrow=N, ncol=nkeep)			# to store posterior beta(t)
	kept_beta_knots 	<- lapply(1:nbeta, matrix, data = NA, nrow=nknots, ncol=nkeep)		# to store posterior beta(t_knots)
	kept_beta_unique 	<- lapply(1:nbeta, matrix, data = NA, nrow=nunique, ncol=nkeep)		# to store posterior beta(t_unique)
	kept_beta_new		<- lapply(1:nbeta, matrix, data = NA, nrow=nnew, ncol=nkeep)		# to store posterior beta(t_new)
	kept_tau 			<- lapply(1:nbeta, matrix, data = NA, nrow=1, ncol=nkeep)			# to store posterior tau
	kept_phi			<- lapply(1:nbeta, matrix, data = NA, nrow=1, ncol=nkeep)			# to store posterior phi
	kept_alpha 			<- matrix(NA,nalpha,nkeep)
	kept_Y 				<- matrix(NA, N, nkeep)
	kept_sigma			<- matrix(NA, 1, nkeep)

	iter 				<- 1																# iteration counter
	ikeep 				<- 1																# kept sample counter
	niter 				<- nburn + nthin*nkeep												# total iterations
																	
	progressbar 		<- progress_bar$new(												# monitor MCMC progress
							format = "[:bar] :current/:total :elapsedfull |:eta",
							total = niter, 												  
							clear = FALSE, 													
							width= 70)

	# progressbar 		<- progressBar(min = 1, max = niter, 
	# 						initial = 1, 
	# 						style = "ETA", 
	# 						substyle = 1,
	# 						char = "#", 
	# 						width = 75, 
	# 						file = "")

	while(iter <= niter){

		#---------------------------- Update Loop ----------------------------#

		# -$-$-$-$-$-$-$-$- update latent Y -$-$-$-$-$-$-$-$-
		if(max(cj)==1){
			h 				<- log(Y_pool)
		}else{
			Y 				<- sampleLatentBM(Y, Y_pool, Xbeta + Galpha, I, Z, sigma, npools, N)
			h 				<- log(Y)
		}
		# -$-$-$-$-$-$-$-$- update beta -$-$-$-$-$-$-$-$-
 		for (d in 1:nbeta){
			Xd				<- X[,d]														# d_th column in X fixed effect desgin matrix
			Xbeta_minus_d 	<- rep(0,N)
			Xbeta_minus_d	<- createXbeta(N,nbeta-1,Xbeta_minus_d,
												matrix(X[,-d],N,nbeta-1),
												matrix(beta[,-d],N,nbeta-1)) 				# Xbeta without dth covariate and dth beta
			h_beta			<- h - Xbeta_minus_d - Galpha									# construct h_beta = h - h without beta_d
			E_mat 			<- gpp_config$E_mat 											# beta_d indicator mat, from t_unique to t seq
			Q_mat 			<- Q_list[[d]] 													# beta_d transform mat, from t_knots to t_unique
			R_knots_mat 	<- R_knots_list[[d]] 											# beta_d knots correlation matrix 
			R_knots_inv_mat <- R_knots_inv_list[[d]]										# beta_d inverse knots correlation matrix
			R_cross_mat 	<- R_cross_list[[d]]											# beta_d nunique x nknots cross corr matrix
			beta_d_config   <- beta_config[,d]												# beta_d config matrix
			
			# update beta_d fun and pour ouput into a list
			beta_d_output 	<- bayesGppFun(center=(d==1), 
											Xd=Xd, sigma=sigma, h=h, h_beta=h_beta, 
											E_mat=E_mat, Q_mat=Q_mat, 
											knots_dist=knots_dist, cross_dist=cross_dist, 
											R_knots_mat=R_knots_mat, 
											R_knots_inv_mat=R_knots_inv_mat, 
											R_cross_mat=R_cross_mat,
											beta_d_config=beta_d_config)

			beta[,d] 				<- beta_d_output$beta_d 								# update beta_d(t)
			beta_unique[,d] 		<- beta_d_output$beta_d_unique							# update beta_d(t_unique)
			beta_knots[,d] 			<- beta_d_output$beta_d_knots 							# update beta_d(t_knots)
			beta_config[3,d] 		<- beta_d_output$tau 									# update tau in beta_config (d_th)
			beta_config[6,d] 		<- beta_d_output$phi 									# update phi in beta_config (d_th)
			beta_config[7,d] 		<- beta_d_output$trphi									# update trphi in beta_config (d_th)
			R_knots_list[[d]] 		<- beta_d_output$R_knots_mat 							# update knots corr mat (d_th)
			R_knots_inv_list[[d]] 	<- beta_d_output$R_knots_inv_mat						# update inverse knots corr mat (d_th)
			R_cross_list[[d]] 		<- beta_d_output$R_cross_mat							# update nunique x nknots cross corr mat (d_th)
			Q_list[[d]] 			<- beta_d_output$Q_mat 									# update transform mat (d_th)
 		
 		}

 		Xbeta 			<- rep(0,N)
		Xbeta 			<- createXbeta(N,nbeta,Xbeta,X,beta)

	    # -$-$-$-$-$-$-$-$- update linear effect alpha -$-$-$-$-$-$-$-$-

 		Sigma_alpha 		<- solve(t(G)%*%G)
 		mu_alpha 			<- Sigma_alpha%*%colSums(G*(h-Xbeta))
 		#alpha				<- as.vector(rmvnorm(1,mu_alpha,as.matrix(Sigma_alpha),method = "svd"))
 		alpha 				<- (t(chol(Sigma_alpha*(sigma**2))))%*%rnorm(nalpha) + mu_alpha
 		Galpha 				<- as.numeric(G%*%alpha)

 		# -$-$-$-$-$-$-$-$- update variance sigma -$-$-$-$-$-$-$-$-

 		sigma 				<- sqrt(1/rgamma(1,N/2+0.5,sum((h-Xbeta-Galpha)^2)/2+0.5))

	    #---------------------------- Save Loop ----------------------------#

	    if(iter>nburn&iter%%nthin==0){

	    	for (d in 1:nbeta){											
	    																# save beta related samples
	    		kept_beta[[d]][,ikeep] 			<- beta[,d]
	    		kept_beta_knots[[d]][,ikeep] 	<- beta_knots[,d]
	    		kept_beta_unique[[d]][,ikeep] 	<- beta_unique[,d]
	    		kept_tau[[d]][,ikeep]			<- beta_config[3,d]	 
	    		kept_phi[[d]][,ikeep]			<- beta_config[6,d]

	    																# save beta(t_new) samples
				
				R_knots_mat_new 				<- matern(cross_dist_new, phi=beta_config[6,d], kappa=beta_config[9,d])
				Q_mat_new 						<- R_knots_mat_new%*%R_knots_inv_list[[d]]
				Sigma_mat 						<- rep(1,nnew)-diag(tcrossprod(Q_mat_new,R_knots_mat_new))
				Sigma_mat[Sigma_mat<=0] 		<- 0
				kept_beta_new[[d]][,ikeep] 		<- as.vector(rnorm(nnew, Q_mat_new%*%beta_knots[,d], Sigma_mat/beta_config[3,d]))
			
			}

																		# save alpha
			kept_alpha[,ikeep]					<- alpha
																		# save latent Y
			kept_Y[,ikeep]						<- Y

																		# save sigma
			kept_sigma[,ikeep]					<- sigma

	    	#---------------------------- Visualization Block ----------------------------#
	    	
	    	# initialize empty moving median storing variables
	    	if(ikeep==1){
					Y11_med 				<- c()
					beta81_med				<- c()
					alpha1_med 				<- c()
					sigma_med 				<- c()
			}
			
			if(visual){
				Y11_med 	<- c(Y11_med,median(kept_Y[1,1:ikeep]))
				beta81_med 	<- c(beta81_med,median(kept_beta[[2]][8,1:ikeep]))
		    	alpha1_med 	<- c(alpha1_med,median(kept_alpha[1,1:ikeep]))
		    	sigma_med 	<- c(sigma_med,median(kept_sigma[1,1:ikeep]))
			}

		    if(visual&iter>nburn&iter%%nthin==0&ikeep%%nstep==0){
		    	
				MC_seq 	<- c(kept_Y[1,1:ikeep],
							kept_beta[[2]][8,1:ikeep],
							kept_alpha[1,1:ikeep],
							kept_sigma[1,1:ikeep])
				
				MC_med <- c(Y11_med, 
							beta81_med,
							alpha1_med, 
							sigma_med)

				MC_true <- c(rep(Y_indv[1],ikeep),
							rep(beta_true[8,2],ikeep),
							rep(alpha_true[1],ikeep),
							rep(sigma_true,ikeep))

				MC_label <- factor(rep(c('Y11',
										'beta81',
										'alpha1',
										'sigma'),each=ikeep))

				MC_df <- data.frame(true=MC_true,values=MC_seq, median=MC_med, iters=rep(1:ikeep,4), label=MC_label)
				
				filename <- paste0('./output/',
									paste0(rep('0',nchar(nkeep/nstep)-nchar(ikeep/nstep)),collapse=''),
									ikeep/nstep,'.png')
				
				trellis.device(png,width=540,height=527,file=filename)
				
				MC_figure <-xyplot(values+true+median~iters|label,
								#main='Markov Chain Monte Carlo',
								xlab=paste0(iter,'th Iteration, ',
											'Draw ',round(ikeep/nkeep*100,2),'% out of ',nkeep,' Samples'),
								ylab='MCMC',
								auto.key=list(space='top',
						 					border=FALSE, 
						 					points=F, 
						 					lines=T, 
						 					columns = 3, 
						 					text=c('MC Chain','True','Moving Median')),
			    				data=MC_df,type="l",
			    				par.settings =list(superpose.line=list(col=c("gray","#28B463",'#A93226'), lwd = c(1.2,3,3), lty=c(1,4,2))),
			    				scales=list(y=list(relation="free")),
			    				strip=strip.custom(factor.levels=expression(alpha[1],
			    															beta[81],
			    															sigma,
			    															Y[11])),
			    				layout=c(2,2))	
				
				print(MC_figure)
				
				dev.off()
				print(MC_figure)
		    }

			ikeep <- ikeep + 1

	    }

	    #setTxtProgressBar(progressbar, iter)
		progressbar$tick()													# set up progress tick
	    
	    iter <- iter + 1

	}

	if(visual){	
		system("convert -delay 50 -loop 0 ./output/*.png ./output/" %+% paste0(type,"mc.gif"))
		system("rm -rf ./output/*.png")
	}

	#------------------------------------- End Gibbs sampling ------------------------------------#

	output_list					<- list(t_knots=t_knots, t_new=t_new, beta_new=kept_beta_new, 
										alpha=kept_alpha, Y=kept_Y, sigma=kept_sigma,
									  	beta_knots=kept_beta_knots, beta_unique=kept_beta_unique, beta=kept_beta,
									  	tau=kept_tau, phi=kept_phi)
		
	return(output_list)

}