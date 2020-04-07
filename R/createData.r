createData <- function(npools=1000, size=5, weighted=TRUE, model=2, copula_rho=0.2,
						default_config=list(u_lower=-2,u_upper=2,
											vol_lower=1,vol_upper=10,
											sigma_true=0.5)){
	cj 			<- rep(size,npools)    
	N 			<- sum(cj)        		
	vol_lower 	<- default_config$vol_lower
	vol_upper 	<- default_config$vol_upper
	u_lower 	<- default_config$u_lower
	u_upper 	<- default_config$u_upper
	sigma_true	<- default_config$sigma_true

	nctn 		<- 2
	copula_mat	<- matrix(copula_rho,nctn,nctn)-diag(copula_rho,nctn)+diag(1,nctn)
	mvnorm_var 	<- rmvnorm(N,sigma=copula_mat)
	X 			<- cbind(1,mvnorm_var[,2:nctn])
	G 			<- cbind(1,rnorm(N),rbinom(N,1,0.4))
	alpha_true 	<- c(-0.5,0.5,0.5)
	unif_var	<- pnorm(mvnorm_var[,1])
	u 			<- unif_var*(u_upper-u_lower)+u_lower # age
	X_homo 		<- X[order(u),]
	G_homo 		<- G[order(u),]
	u_homo 		<- u[order(u)]
	W_mat 		<- matrix(NA,npools,max(cj)) # weight matrix
	vol   		<- runif(N, vol_lower, vol_upper) # volume for each specimen
	vol_homo 	<- vol[order(u)]

	logit 		<<- function(x) log(x/(1-x))
	logit_inv 	<<- function(x) return(1/(1+exp(-x)))


	if((u_lower==-2)&(u_upper=2)){
	  # Define beta functions, age range [-2,2]
	  f1 <- function(u,k=1.2,sig=2,a=2){
	    0.7*exp((-(u>0)*k^a-(u<0)/k^a)*abs(u/sig)^a)
	  }
	  
	  f2 <- function(u,mu1=-1,mu2=1,sig1=.6,sig2=.7,wt=.6){
	    wt*exp(-(u-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(u-mu2)^2/(2*sig2^2))
	  }
	  
	  f3<- function(u){
	    0.2*(sin(pi*(u-0.5)/2.5)+2)*(exp(-(u+(u+0.5)^2*(u>0.5))/6))
	  }
	  
	  f4 <- function(u,a=0,b=2,c=4){
	    c*exp(a+b*u)/(6+6*exp(a+b*u))
	  }
	  
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
	}else{
	  # Define beta functions, age range [-3,3]
	  f1 <- function(u,k=1.2,sig=2.5,a=2){
	    0.7*exp((-(u>0)*k^a-(u<0)/k^a)*abs(u/sig)^a)
	  }
	  
	  f2 <- function(u,mu1=-1.5,mu2=1.5,sig1=.6,sig2=.8,wt=.6){
	    wt*exp(-(u-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(u-mu2)^2/(2*sig2^2))
	  }
	  
	  f3<- function(u){
	    0.2*(sin(pi*(u+0.3)/2.5)+2)*(exp(-(u+(u-0.3)^2*(u>0.3))/6))
	  }
	  
	  f4 <- function(u,a=1,b=1.5,c=4){
	    c*exp(a+b*u)/(6+6*exp(a+b*u))
	  }
	  
	  betaTrue <<- function(u, beta_list, center_intc=TRUE){
	    u_lower 		<- min(u); u_upper <- max(u)
	    N 				<- length(u)
	    nbeta 			<- length(beta_list) # total number of beta: beta0 beta1 beta2 beta3
	    beta 			<- matrix(NA, N, nbeta)
	    
	    for (d in 1:nbeta){
	      constant 		<- integrate(beta_list[[d]],u_lower,u_upper)$value/(u_upper-u_lower)
	      beta[,d] 		<- beta_list[[d]](u) - constant*(d==(1*center_intc))
	      #beta[,d] 		<- beta_list[[d]](u) - constant*center_intc
	    }
	    
	    return(beta)
	  }
	}

	if(model == 1){ # combine of function 1 and 2
	  
	  beta_list 	<- list(beta_0=f1, beta_1=f2)
	  beta_true 	<- betaTrue(u=u,beta_list=beta_list,center_intc=TRUE)
	  beta_homo_true<- beta_true[order(u),]
	  Y_indv 		<- rlnorm(N,G%*%alpha_true+apply(X*beta_true,1,sum),sdlog=sigma_true)
	  Y_homo_indv 	<- Y_indv[order(u)]
	  
	}else if(model == 2){ # combine of function 1 and 2
	  
	  beta_list 	<- list(beta_0=f3, beta_1=f4)
	  beta_true 	<- betaTrue(u=u,beta_list=beta_list,center_intc=TRUE)
	  beta_homo_true<- beta_true[order(u),]
	  Y_indv 		<- rlnorm(N,G%*%alpha_true+apply(X*beta_true,1,sum),sdlog=sigma_true)
	  Y_homo_indv 	<- Y_indv[order(u)]
	  
	}

	if(weighted){
	  vol_split 	<- split(vol,rep(1:npools,cj))
	  wgt_lst 		<- lapply(vol_split,function(vec) vec/sum(vec))
	  wgt_vec 		<- unlist(wgt_lst)
	  W_mat 		<- matrix(wgt_vec,nrow=npools,byrow=T)
	  
	  vol_homo_split<- split(vol_homo,rep(1:npools,cj))
	  wgt_homo_lst 	<- lapply(vol_homo_split,function(vec) vec/sum(vec))
	  wgt_homo_vec 	<- unlist(wgt_homo_lst)
	  W_homo_mat 	<- matrix(wgt_homo_vec,nrow=npools,byrow=T)
	}else{
	  vol_split 	<- split(rep(1,N),rep(1:npools,cj))
	  wgt_lst 		<- lapply(vol_split,function(vec) vec/sum(vec))
	  wgt_vec 		<- unlist(wgt_lst)
	  W_mat 		<- matrix(wgt_vec,nrow=npools,byrow=T)
	  
	  vol_homo_split<- split(rep(1,N),rep(1:npools,cj))
	  wgt_homo_lst 	<- lapply(vol_homo_split,function(vec) vec/sum(vec))
	  wgt_homo_vec 	<- unlist(wgt_homo_lst)
	  W_homo_mat 	<- matrix(wgt_homo_vec,nrow=npools,byrow=T)
	}

	pool_id 		<- rep(1:npools, cj) # latent invidual-wise
	I_mat 			<- cbind(0,pool_id,as.vector(t(W_mat)))
	I_homo_mat 		<- cbind(0,pool_id,as.vector(t(W_homo_mat)))

	Y_pool      	<- aggregate(Y_indv*I_mat[,3]~I_mat[,2],FUN=sum)[,2] # observed pool-wise
	Y_homo_pool 	<- aggregate(Y_homo_indv*I_homo_mat[,3]~I_homo_mat[,2],FUN=sum)[,2]

	Z_mat    		<- matrix(NA,npools,2+max(cj))
	Z_homo_mat 		<- matrix(NA,npools,2+max(cj))
	ind_1st 		<- 0

	for(j in 1:npools){
		Z_mat[j,]   	<- c(Y_pool[j],cj[j],(ind_1st+1):(ind_1st+cj[j]))
		Z_homo_mat[j,]	<- c(Y_homo_pool[j],cj[j],(ind_1st+1):(ind_1st+cj[j]))
		ind_1st 			<- ind_1st+cj[j]
	}

	output_lst <- list(Y_indv = Y_indv, Y_pool=Y_pool, I_mat=I_mat, Z_mat=Z_mat,
						W_mat=W_mat, X=X, u=u, G=G, alpha_true=alpha_true, beta_true=beta_true, sigma_true=sigma_true,
						Y_homo_indv = Y_indv, Y_homo_pool=Y_homo_pool, I_homo_mat=I_mat, Z_homo_mat=Z_mat, 
						W_homo_mat=W_mat, X_homo=X, u_homo=u_homo, G_homo=G_homo, beta_homo_true=beta_homo_true)

	return(output_lst)
}