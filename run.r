rm(list=ls(all=TRUE))
graphics.off()
#options(echo=FALSE)
options(echo=TRUE)
cluster <- FALSE
keys <- c('taskid','reps','nburn','nkeep','nthin','nstep','npools','size','weighted','copula_rho','model','homo','plot','name')
for(key in keys) assign(key, NA)
input_args <- commandArgs(trailingOnly = TRUE)

if(!cluster){
	setwd("/Users/Zeng/Dropbox/research-at-sc/vcm/vcm-biomarker")
	input_args <- c('taskid=1',
					'reps=1',
					'nburn=1000',
					'nkeep=2000',
					'nthin=10',
					'nstep=10',
					'npools=500',
					'size=5',
					'weighted=TRUE',
					'copula_rho=0.2',
					'model=1',
					'homo=FALSE',
					'plot=TRUE',
					"name='biomaker'")
}

if(length(input_args)!=0){
	input_dict = unlist(strsplit(input_args,'='))
	input_keys = input_dict[seq(1,length(input_dict),2)]
	if(!all(input_keys %in% keys)) stop('[error] incorrect arguments names!')
	eval(parse(text=input_args))
}

if(identical(taskid, NA)) taskid <- 1
if(identical(reps, NA)) reps <- 1
if(identical(nburn, NA)) nburn <- 1000
if(identical(nkeep, NA)) nkeep <- 2000
if(identical(nthin, NA)) nthin <- 5
if(identical(nstep, NA)) nstep <- 20
if(identical(npools, NA)) plot <- 500
if(identical(size, NA)) size <- 5
if(identical(weighted, NA)) weighted <- TRUE
if(identical(copula_rho, NA)) copula_rho <- 0.2
if(identical(model, NA)) model <- 1
if(identical(homo, NA)) homo <- FALSE
if(identical(plot, NA)) plot <- FALSE
if(identical(name, NA)) name <- 'biomaker'

source('./R/libraries.r')
for(rep in 1:reps){
	
	#output_file = './output/output' %+% as.character(reps*(taskid-1)+rep) %+% '.md'
	#if (file.exists(output_file)) file_remove = file.remove(output_file)
	#cat(paste0(paste0(rep('=',10),collapse = ''),'Bayesian VCM with BM',paste0(rep('=',10),collapse = '')),file=output_file)
	
	data_lst 	<- createData(npools=npools, size=size, weighted=weighted, model=model, copula_rho=copula_rho)

	X 			<- data_lst$X
	G 			<- data_lst$G
	u 			<- data_lst$u
	Y_obs 		<- data_lst$Y_pool
	I_mat 		<- data_lst$I_mat
	Z_mat 		<- data_lst$Z_mat
	W_mat 		<- data_lst$W_mat

	X_homo 		<- data_lst$X_homo
	G_homo		<- data_lst$G_homo
	u_homo 		<- data_lst$u_homo
	Y_homo_obs 	<- data_lst$Y_homo_pool
	I_homo_mat 	<- data_lst$I_homo_mat
	Z_homo_mat 	<- data_lst$Z_homo_mat
	W_homo_mat 	<- data_lst$W_homo_mat

	Y_indv 		<- data_lst$Y_indv
	Y_homo_indv <- data_lst$Y_homo_indv
	alpha_true 	<- data_lst$alpha_true
	beta_true 	<- data_lst$beta_true
	beta_homo_true 	<- data_lst$beta_homo_true
	sigma_true 	<- data_lst$sigma_true

	true_params			<- list(Y_indv=Y_indv,sigma_true=sigma_true,beta_true=beta_true,alpha_true=alpha_true)
	true_homo_params	<- list(Y_indv=Y_homo_indv,sigma_true=sigma_true,beta_true=beta_homo_true)
	
	# if(reps*(taskid-1)+rep==1){
	# 	postfix <- 'st'
	# }else if(reps*(taskid-1)+rep==2){
	# 	postfix <- 'nd'
	# }else if(reps*(taskid-1)+rep==3){
	# 	postfix <- 'rd'
	# }else{
	# 	postfix <- 'th'
	# }
	#header <- paste0('\n','[Done!] The ',paste0(reps*(taskid-1)+rep,postfix)) %+% paste0(' simulation with prevalence ',format(round(mean(Y_true),digits=2),nsmall=2),', based on')

	out_BM 			<- try(bayesGppBM(X=X, G=G, t=u, Y_pool=Y_obs, Z=Z_mat, I=I_mat, W=W_mat, size=size, 
							visual=plot, true_params=true_params,
							hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
							gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)

	if(homo){
		out_BM 			<- try(bayesGppBM(X=X, G=G, t=u, Y_pool=Y_obs, Z=Z_mat, I=I_mat, W=W_mat, size=size, 
							visual=plot, true_params=true_homo_params,
							hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
							gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
	}
	#saveData(MC_list=out_MPT, taskid=taskid, reps=reps, rep=rep, savepath='./output/MPT/')
	#MPT_cue 					<- header %+% ', MPT'
	#cat(MPT_cue,file=output_file,append=TRUE)
	
}


























