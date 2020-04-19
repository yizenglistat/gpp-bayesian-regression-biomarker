createData <- function(npools=500, sizes=c(1,2,4,6,8,10), weighted=TRUE, model=2, copula_rho=0.2,
						default_config=list(u_lower=-2,u_upper=2,
											vol_lower=1,vol_upper=10,
											sigma_true=0.5)){

	if(TRUE){
		npools=500
		sizes=c(1,2,4,6,8,10)
		weighted=TRUE;
		model=2;
		copula_rho=0.2;
		default_config=list(u_lower=-2,u_upper=2,vol_lower=1,vol_upper=10,sigma_true=0.5)
	}

	logit 		<<- function(x) log(x/(1-x))
	logit_inv 	<<- function(x) return(1/(1+exp(-x)))
	
	betaTrue <<- function(u, beta_list, center_intc=TRUE){
		u_lower 		<- min(u);	u_upper <- max(u)
		N 				<- length(u)
		nbeta 			<- length(beta_list) # total number of beta: beta0 beta1 beta2 beta3
		beta 			<- matrix(NA, N, nbeta)

		for (d in 1:nbeta){
		  const 		<- integrate(beta_list[[d]],u_lower,u_upper)$value/(u_upper-u_lower)
		  # depd_const	<- integrate(function(u) beta_list[[d]](u)*dnorm(u),-15,15)$value
		  beta[,d] 	<- beta_list[[d]](u) - const*(d==(1*center_intc))
		}

		return(beta)
	}

	model_list				<-list()

 	# Define beta functions, age range [-2,2]
	f1 <- function(u, scale=1/8) scale*u^3
	f2 <- function(u, shape=pi/2) sin(shape*u)
	f3 <- function(u, scale=0.25) scale*u*(1-u)
	f4 <- function(u, shape=-1) exp(shape*u^2) 
	#model_list[['model1']] <- list(beta_0=f1, beta_1=f2)
	model_list[['model1']] <- list(beta_0=f1, beta_1=f2, beta_2=f3, beta_3=f4)
	
	# Define beta functions, age range [-2,2]
	f1 <- function(u, scale=1/2) scale*u^2
	f2 <- function(u, shape=1.44,scale=4) exp((-1*I(u>0)*shape-I(u<0)/shape)*u^2/scale)
	f3 <- function(u, rate1=0.72,rate2=0.98,scale=0.7) exp(-(u+1)^2/rate1)+scale*exp(-(u-1)^2/rate2)
	f4 <- function(u, a=0.4,b=0.8,c=2.5,d=6) (a*sin(pi*(u-0.5)/c)+b)/exp((u+(u+0.5)^2*I(u>-0.5))/d) 
	#model_list[['model2']] <- list(beta_0=f1, beta_1=f2)
	model_list[['model2']] <- list(beta_0=f1, beta_1=f2, beta_2=f3, beta_3=f4)

	cj 			<- rep(max(sizes),npools)    
	N 			<- sum(cj)        		
	vol_lower 	<- default_config$vol_lower
	vol_upper 	<- default_config$vol_upper
	u_lower 	<- default_config$u_lower
	u_upper 	<- default_config$u_upper
	sigma_true	<- default_config$sigma_true

	nctn 		<- 2
	copula_mat	<- matrix(copula_rho,nctn,nctn)-diag(copula_rho,nctn)+diag(1,nctn)
	mvnorm_var 	<- rmvnorm(N,sigma=copula_mat)
	X 			<- cbind(1,mvnorm_var[,2:nctn],rbinom(N,1,0.5),rbinom(N,1,0.3))
	unif_var	<- pnorm(mvnorm_var[,1])
	u 			<- round(unif_var*(u_upper-u_lower)+u_lower,2) # age
	vol   		<- runif(N, vol_lower, vol_upper) # volume for each specimen

	if(model == 1){ # combine of function 1 and 2
		beta_list 	<- model_list$model1
		beta_true 	<- betaTrue(u=u,beta_list=beta_list,center_intc=FALSE)
		Y_indv 		<- rlnorm(N, apply(X*beta_true,1,sum), sdlog=sigma_true)
	  	  
	}else if(model == 2){ # combine of function 1 and 2

	 	beta_list 	<- model_list$model2
	  	beta_true 	<- betaTrue(u=u,beta_list=beta_list,center_intc=FALSE)
	 	Y_indv 		<- rlnorm(N,apply(X*beta_true,1,sum),sdlog=sigma_true)
	  
	}
	
	sortData <- function(X,u,vol,beta_true,Y_indv,size,ngroups=48){
		
		# homo sort
		homo_order 		<- order(X[,3],X[,4],u) # order by age, gender, race non-decreasing	
		X_homo 			<- X[homo_order,]
		u_homo 			<- u[homo_order]
		vol_homo 		<- vol[homo_order]
		beta_homo_true 	<- beta_true[homo_order,]
		Y_homo_indv 	<- Y_indv[homo_order]

		# nearly homo sort
		cj 				<- rep(size,npools) 
		N 				<- nrow(X)
		idx_split 		<- list()
		for(group in 1:ngroups){
			group_size <- floor(N/ngroups)
			idx_split[[group]] <- (1+(group-1)*group_size):( (group*group_size)*(group!=48)+ N*(group==48))
		}
		idx_lst 		<- lapply(idx_split,function(idx) idx[order(vol_homo[idx])])
		nhomo_order 	<- unlist(idx_lst)
		X_nhomo 		<- X_homo[nhomo_order,]
		u_nhomo 		<- u_homo[nhomo_order]
		vol_nhomo 		<- vol_homo[nhomo_order]
		beta_nhomo_true <- beta_homo_true[nhomo_order,]
		Y_nhomo_indv 	<- Y_homo_indv[nhomo_order]

		if(weighted){

			vol_split 		<- split(vol,rep(1:npools,cj))
			wgt_lst 		<- lapply(vol_split,function(vec) vec/sum(vec))
			wgt_vec 		<- unlist(wgt_lst)
			W_mat 			<- matrix(wgt_vec,nrow=npools,byrow=TRUE)

			vol_homo_split 	<- split(vol_homo,rep(1:npools,cj))
			wgt_homo_lst 	<- lapply(vol_homo_split,function(vec) vec/sum(vec))
			wgt_homo_vec 	<- unlist(wgt_homo_lst)
			W_homo_mat 		<- matrix(wgt_homo_vec,nrow=npools,byrow=TRUE)

			vol_nhomo_split <- split(vol_nhomo,rep(1:npools,cj))
			wgt_nhomo_lst 	<- lapply(vol_nhomo_split,function(vec) vec/sum(vec))
			wgt_nhomo_vec 	<- unlist(wgt_nhomo_lst)
			W_nhomo_mat 	<- matrix(wgt_nhomo_vec,nrow=npools,byrow=TRUE)

		}else{

			vol_split 		<- split(rep(1,N),rep(1:npools,cj))
			wgt_lst 		<- lapply(vol_split,function(vec) vec/sum(vec))
			wgt_vec 		<- unlist(wgt_lst)
			W_mat 			<- matrix(wgt_vec,nrow=npools,byrow=TRUE)

			vol_homo_split 	<- split(rep(1,N),rep(1:npools,cj))
			wgt_homo_lst 	<- lapply(vol_homo_split,function(vec) vec/sum(vec))
			wgt_homo_vec 	<- unlist(wgt_homo_lst)
			W_homo_mat 		<- matrix(wgt_homo_vec,nrow=npools,byrow=TRUE)

			vol_nhomo_split <- split(rep(1,N),rep(1:npools,cj))
			wgt_nhomo_lst 	<- lapply(vol_nhomo_split,function(vec) vec/sum(vec))
			wgt_nhomo_vec 	<- unlist(wgt_nhomo_lst)
			W_nhomo_mat 	<- matrix(wgt_nhomo_vec,nrow=npools,byrow=TRUE)

		}

		pool_id 		<- rep(1:npools, cj) # latent invidual-wise
		I_mat 			<- cbind(0,pool_id,as.vector(t(W_mat)))
		I_homo_mat 		<- cbind(0,pool_id,as.vector(t(W_homo_mat)))
		I_nhomo_mat 	<- cbind(0,pool_id,as.vector(t(W_nhomo_mat)))

		Y_pool      	<- aggregate(Y_indv*I_mat[,3]~I_mat[,2],FUN=sum)[,2] # observed pool-wise
		Y_homo_pool 	<- aggregate(Y_homo_indv*I_homo_mat[,3]~I_homo_mat[,2],FUN=sum)[,2]
		Y_nhomo_pool 	<- aggregate(Y_nhomo_indv*I_nhomo_mat[,3]~I_nhomo_mat[,2],FUN=sum)[,2]

		Z_mat    		<- matrix(NA,npools,2+max(cj))
		Z_homo_mat 		<- matrix(NA,npools,2+max(cj))
		Z_nhomo_mat 	<- matrix(NA,npools,2+max(cj))
		ind_1st 		<- 0

		for(j in 1:npools){
			Z_mat[j,]   	<- c(Y_pool[j],cj[j],(ind_1st+1):(ind_1st+cj[j]))
			Z_homo_mat[j,]	<- c(Y_homo_pool[j],cj[j],(ind_1st+1):(ind_1st+cj[j]))
			Z_nhomo_mat[j,]	<- c(Y_nhomo_pool[j],cj[j],(ind_1st+1):(ind_1st+cj[j]))
			ind_1st 		<- ind_1st+cj[j]
		}

		output_lst <- list(Y_indv = Y_indv, Y_pool=Y_pool, I_mat=I_mat, Z_mat=Z_mat, model_list=model_list,
							W_mat=W_mat, X=X, u=u, beta_true=beta_true, sigma_true=sigma_true,
							
							Y_homo_indv = Y_homo_indv, Y_homo_pool=Y_homo_pool, I_homo_mat=I_homo_mat, Z_homo_mat=Z_homo_mat, 
							W_homo_mat=W_homo_mat, X_homo=X_homo, u_homo=u_homo, beta_homo_true=beta_homo_true,

							Y_nhomo_indv = Y_nhomo_indv, Y_nhomo_pool=Y_nhomo_pool, I_nhomo_mat=I_nhomo_mat, Z_nhomo_mat=Z_nhomo_mat, 
							W_nhomo_mat=W_nhomo_mat, X_nhomo=X_nhomo, u_nhomo=u_nhomo, beta_nhomo_true=beta_nhomo_true)
		return(output_lst)
	}

	output_lst <- list()
	for(size in sizes){

		idx_size 		<- 1:(npools*size)
		X_size 			<- X[idx_size,]
		u_size 			<- u[idx_size]
		vol_size 		<- vol[idx_size]
		beta_true_size 	<- beta_true[idx_size,]
		Y_indv_size 	<- Y_indv[idx_size]
		
		output_lst[[size]] <- sortData(X=X_size,
									u=u_size,
									vol=vol_size,
									beta_true=beta_true_size,
									Y_indv=Y_indv_size,size=size)
		
	}

	return(output_lst)

}