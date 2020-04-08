rm(list=ls(all=TRUE))
graphics.off()
#options(echo=FALSE)
options(echo=TRUE)
cluster <- FALSE
keys <- c('taskid','reps','nburn','nkeep','nthin','nstep','npools','sizes','weighted','copula_rho','models','homo','plot','name')
for(key in keys) assign(key, NA)
input_args <- commandArgs(trailingOnly = TRUE)

if(!cluster){
	setwd("/Users/Zeng/Dropbox/research-at-sc/vcm/vcm-biomarker")
	input_args <- c('taskid=1',
					'reps=1',
					'nburn=1',
					'nkeep=5',
					'nthin=2',
					'nstep=2',
					'npools=300',
					'sizes=c(1,2,4,6,8)',
					'weighted=TRUE',
					'copula_rho=0.2',
					'models=c(1,2)',
					'homo=TRUE',
					'plot=FALSE',
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
if(identical(nburn, NA)) nburn <- 5000
if(identical(nkeep, NA)) nkeep <- 2000
if(identical(nthin, NA)) nthin <- 10
if(identical(nstep, NA)) nstep <- 20
if(identical(npools, NA)) plot <- 500
if(identical(sizes, NA)) sizes <- 2
if(identical(weighted, NA)) weighted <- TRUE
if(identical(copula_rho, NA)) copula_rho <- 0.2
if(identical(models, NA)) models <- 1
if(identical(homo, NA)) homo <- FALSE
if(identical(plot, NA)) plot <- FALSE
if(identical(name, NA)) name <- 'biomaker'

source('./R/libraries.r')
output_file = './output/output.md'
if (file.exists(output_file)) file_remove = file.remove(output_file)
cat(paste0(paste0(rep('=',30),collapse = ''),'Bayesian VCM with BM',paste0(rep('=',30),collapse = '')),file=output_file,'\n')

for(model in models){

	for(size in sizes){
		
		for(rep in 1:reps){

			if(reps*(taskid-1)+rep==1){
				postfix <- 'st'
			}else if(reps*(taskid-1)+rep==2){
				postfix <- 'nd'
			}else if(reps*(taskid-1)+rep==3){
				postfix <- 'rd'
			}else{
				postfix <- 'th'
			}

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

			Y_indv 		<- data_lst$Y_indv
			Y_homo_indv <- data_lst$Y_homo_indv
			alpha_true 	<- data_lst$alpha_true
			beta_true 	<- data_lst$beta_true
			beta_homo_true 	<- data_lst$beta_homo_true
			sigma_true 	<- data_lst$sigma_true
			model_list 	<- data_lst$model_list

			true_params			<- list(Y_indv=Y_indv,sigma_true=sigma_true,beta_true=beta_true,alpha_true=alpha_true)
			true_homo_params	<- list(Y_indv=Y_homo_indv,sigma_true=sigma_true,beta_true=beta_homo_true,alpha_true=alpha_true)
			
			# ========================================= Bayesian Gpp Fitting ========================================= # 
		
			out_BM 				<- try(bayesGppBM(X=X, G=G, t=u, Y_pool=Y_obs, Z=Z_mat, I=I_mat, size=size, 
									visual=plot, true_params=true_params, model=model,
									hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
									gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
			
			saveData(MC_list=out_BM, taskid=taskid, reps=reps, rep=rep, savepath=paste0('./output/model',model,'/size',size,'/'))
		
			if(homo){
				if(size==1){
					out_homo_BM <- out_BM
				}else{
					out_homo_BM 	<- try(bayesGppBM(X=X_homo, G=G_homo, t=u_homo, Y_pool=Y_homo_obs, Z=Z_homo_mat, I=I_homo_mat, size=size, 
									visual=plot, true_params=true_homo_params, model=model,
									hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
									gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
				}

				saveData(MC_list=out_homo_BM, taskid=taskid, reps=reps, rep=rep, savepath=paste0('./output/model',model,'/size',size,'/homo/'))
			}

			header <- paste0('\n','[Done!] The ',paste0(reps*(taskid-1)+rep,postfix)) %+% paste0(' simulation with model=',model)
			cat(paste0(header,', size=',size),file=output_file,append=TRUE)
		}

	}

}

#bayesPlot(plot_beta=FALSE,plot_sigma=TRUE,which_beta=0,model=1,sizes=c(1))












