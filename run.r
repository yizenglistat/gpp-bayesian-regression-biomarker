rm(list=ls(all=TRUE))
graphics.off()
options(echo=FALSE)
options(echo=TRUE)
cluster <- TRUE
cluster <- FALSE
keys <- c('taskid','reps','nburn','nkeep','nthin','nstep','npools','sizes','weighted','copula_rho','models','random','homo','nhomo','plot','name')
for(key in keys) assign(key, NA)
input_args <- commandArgs(trailingOnly = TRUE)

if(!cluster){
	setwd("/Users/Zeng/Dropbox/research-at-sc/vcm/vcm-continuous")
	input_args <- c('taskid=1',
					'reps=1',
					'nburn=1000',
					'nkeep=2000',
					'nthin=10',
					'nstep=20',
					'npools=300',
					'sizes=c(1,2,4,6,8,10)',
					'weighted=TRUE',
					'copula_rho=0.2',
					'models=1',
					'random=FALSE',
					'homo=FALSE',
					'nhomo=TRUE',
					'plot=TRUE',
					"name='Biomaker'")
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
if(identical(random, NA)) random <- TRUE
if(identical(homo, NA)) homo <- FALSE
if(identical(nhomo, NA)) nhomo <- FALSE
if(identical(plot, NA)) plot <- FALSE
if(identical(name, NA)) name <- 'biomaker'

source('./R/libraries.r')
update_files <- function(output_file = './output/output.md'){
	cat(paste0(paste0(rep('=',12),collapse = ''),'Bayesian VCM with BM',paste0(rep('=',12),collapse = '')),file=output_file,'\n')
	cat('model','|','size','|','type','|','# of theta','|','# of beta',file=output_file,append=TRUE,'\n')
	for(model in c(1,2)){
		for(size in c(1,2,4,6,8,10)){
			outpath = paste0('./output/model',model,'/size',size,'/')
			outpath_homo = paste0('./output/model',model,'/size',size,'/','homo/')
			outpath_nhomo = paste0('./output/model',model,'/size',size,'/','nhomo/')
			count_theta = length(list.files(outpath,pattern='theta_summ[0-9]*.csv'))
			count_beta = length(list.files(outpath,pattern='beta_knots_summ[0-9]*.csv'))
			count_theta_homo = length(list.files(outpath_homo,pattern='theta_summ[0-9]*.csv'))
			count_beta_homo = length(list.files(outpath_homo,pattern='beta_knots_summ[0-9]*.csv'))
			count_theta_nhomo = length(list.files(outpath_nhomo,pattern='theta_summ[0-9]*.csv'))
			count_beta_nhomo = length(list.files(outpath_nhomo,pattern='beta_knots_summ[0-9]*.csv'))
			cat('  ',model,'    ',size,paste0(rep(' ',2-nchar(size)),collapse=''),'  ','random     ',count_theta,paste0(rep(' ',10-nchar(count_theta)),collapse=''), count_beta,file=output_file,append=TRUE,'\n')
			cat('  ',model,'    ',size,paste0(rep(' ',2-nchar(size)),collapse=''),'  ','homo       ',count_theta_homo,paste0(rep(' ',10-nchar(count_theta_homo)),collapse=''), count_beta_homo,file=output_file,append=TRUE,'\n')
			cat('  ',model,'    ',size,paste0(rep(' ',2-nchar(size)),collapse=''),'  ','nhomo      ',count_theta_nhomo,paste0(rep(' ',10-nchar(count_theta_nhomo)),collapse=''), count_beta_nhomo,file=output_file,append=TRUE,'\n')
		}
	}
}



update_files()
for(model in models){

	for(rep in 1:reps){
		
		data_lst <- createData(npools=npools, 
						sizes=sizes, 
						weighted=weighted, 
						model=model, 
						copula_rho=copula_rho)

		for(size in rev(sizes)){
			data_lst_of_size <- data_lst[[size]]
			X 			<- data_lst_of_size$X
			u 			<- data_lst_of_size$u
			Y_obs 		<- data_lst_of_size$Y_pool
			I_mat 		<- data_lst_of_size$I_mat
			Z_mat 		<- data_lst_of_size$Z_mat
			W_mat 		<- data_lst_of_size$W_mat

			X_homo 		<- data_lst_of_size$X_homo
			u_homo 		<- data_lst_of_size$u_homo
			Y_homo_obs 	<- data_lst_of_size$Y_homo_pool
			I_homo_mat 	<- data_lst_of_size$I_homo_mat
			Z_homo_mat 	<- data_lst_of_size$Z_homo_mat

			X_nhomo 		<- data_lst_of_size$X_nhomo
			u_nhomo 		<- data_lst_of_size$u_nhomo
			Y_nhomo_obs 	<- data_lst_of_size$Y_nhomo_pool
			I_nhomo_mat 	<- data_lst_of_size$I_nhomo_mat
			Z_nhomo_mat 	<- data_lst_of_size$Z_nhomo_mat

			Y_indv 			<- data_lst_of_size$Y_indv
			Y_homo_indv 	<- data_lst_of_size$Y_homo_indv
			Y_nhomo_indv 	<- data_lst_of_size$Y_nhomo_indv
			beta_true 		<- data_lst_of_size$beta_true
			beta_homo_true 	<- data_lst_of_size$beta_homo_true
			beta_nhomo_true <- data_lst_of_size$beta_nhomo_true
			sigma_true 		<- data_lst_of_size$sigma_true
			model_list 		<- data_lst_of_size$model_list

			true_params			<- list(Y_indv=Y_indv,sigma_true=sigma_true,beta_true=beta_true)
			true_homo_params	<- list(Y_indv=Y_homo_indv,sigma_true=sigma_true,beta_true=beta_homo_true)
			true_nhomo_params	<- list(Y_indv=Y_nhomo_indv,sigma_true=sigma_true,beta_true=beta_nhomo_true)
			
			# ========================================= Bayesian Gpp Fitting ========================================= # 
			
			# ------------ Random Pooling ------------#
			if(random){
				out_BM 				<- try(bayesGppBM(X=X, G=G, t=u, Y_pool=Y_obs, Z=Z_mat, I=I_mat, size=size, 
										visual=plot, true_params=true_params, model=model,
										hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
										gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
				
				savepath = paste0('./output/model',model,'/size',size,'/')
				#ignored = dir.create(savepath, showWarnings = FALSE)
				saveData(MC_list=out_BM, taskid=taskid, reps=reps, rep=rep, savepath=savepath)
				update_files()
			}
			# ------------ Homogeneous Pooling ------------#
			if(homo){
				if(size==1){
					out_homo_BM <- out_BM
				}else{
					out_homo_BM 	<- try(bayesGppBM(X=X_homo, t=u_homo, Y_pool=Y_homo_obs, Z=Z_homo_mat, I=I_homo_mat, size=size, 
									visual=plot, true_params=true_homo_params, model=model,
									hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
									gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
				}

				savepath = paste0('./output/model',model,'/size',size,'/homo/')
				#ignored = dir.create(savepath, showWarnings = FALSE)
				saveData(MC_list=out_homo_BM, taskid=taskid, reps=reps, rep=rep, savepath=savepath)
			}
			update_files()
			# ------------ nearly Homogeneous Pooling ------------#
			if(nhomo){
				if(size==1){
					out_nhomo_BM <- out_BM
				}else{
					out_nhomo_BM 	<- try(bayesGppBM(X=X_nhomo, t=u_nhomo, Y_pool=Y_nhomo_obs, Z=Z_nhomo_mat, I=I_nhomo_mat, size=size, 
									visual=plot, true_params=true_nhomo_params, model=model,
									hyper_params=list(phi_sd=0.1, kappa=2, nknots=50,sigma_init=1),
									gibbs_config=list(nburn=nburn,nkeep=nkeep,nthin=nthin,nstep=nstep)),silent=TRUE)
				}

				savepath = paste0('./output/model',model,'/size',size,'/nhomo/')
				#ignored = dir.create(savepath, showWarnings = FALSE)
				saveData(MC_list=out_nhomo_BM, taskid=taskid, reps=reps, rep=rep, savepath=savepath)
			}
			update_files()

		}

	}

}

#bayesPlot(plot_beta=TRUE,plot_sigma=TRUE,plot_ise=TRUE,which_beta=0,model=1,sizes=c(2,4))




