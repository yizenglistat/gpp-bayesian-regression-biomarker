
if(TRUE){
  weight=TRUE;J=300;N=sum(cj);L_1=50;L_2=50;lbd=c(-2,-2);rbd=c(2,2)
	beta=c(.5,-.5,.5);sig=0.5;esig=0.1;mod=1
}

# Uniform,Normal,Binary(race),Binary(Gender)  # mod 1
# Uniform,Normal,Binary(race),linear normal   # mod 2

# cj is a vector of length J containing pool sizes
data_sim <- function(weight=TRUE,mod,beta,sig,esig,N,cj,J,lbd,rbd,L_1,L_2,copula.rho){
  
  	lbd_1 <- lbd[1]
	lbd_2 <- lbd[2]
	rbd_1 <- rbd[1]
	rbd_2 <- rbd[2]

  # Generate 2 continuous covariates that are nonlinearly related to E(Y)
	x_1.tem <- rnorm(N)
	x_2.tem <- rnorm(N)
	matrix.tem 	<- matrix(c(1,copula.rho,copula.rho,1),2,2)
	solve.tem 	<- eigen(matrix.tem)
	copula.design<-solve.tem$vectors%*%diag(sqrt(solve.tem$values))
	xxxx<-copula.design%*%rbind(x_1.tem,x_2.tem)
	x_1<-xxxx[1,]
	x_2<-xxxx[2,]
	x_1<-pnorm(x_1)*(rbd_1-lbd_1)+lbd_1
	#cor(x_1,x_2)
	#mean(pnorm(x_1));var(pnorm(x_1))
	#mean(x_1);var(x_1)
	#plot(sort(x_1),1:N)
	#plot(sort(x_2),1:N)
	#plot(x_1,x_2)
	
	#x_1   <- runif(N,lbd_1,rbd_1) # Age
	#x_2   <- rnorm(N)             # BMI
	Wmat  <- cbind(1,rbinom(N,1,0.5),rbinom(N,1,0.2)) # Design matrix
	Omat  <- matrix(NA,J,max(cj))                     # Weight matrix
	Vol   <- runif(N,1,10)                            # Generate volume of each specimen
	
	# ====================================== Homogeneous pool
	Omat_HG  <- matrix(NA,J,max(cj))                     # Weight matrix of Homegeneous pooling
	oid      <- order(Wmat[,2],Wmat[,3],x_1)
	Wmat_HG  <- Wmat[oid,]
	x_1_HG   <- x_1[oid]
	x_2_HG   <- x_2[oid]
	Vol_HG   <- Vol[oid]
	
	#head(cbind(Wmat,x_1),20)
	#head(cbind(Wmat_HG,x_1_HG),20)
	
	csm <- 0
	for(jj in 1:J){
		if(weight){
			x    <- Vol[(csm+1):(csm+cj[jj])]
			x_HG <- Vol_HG[(csm+1):(csm+cj[jj])]
		}else{
			x    <- rep(1,cj[jj])
			x_HG <- rep(1,cj[jj])
		}
		Omat[jj,]    <- x/sum(x)
		Omat_HG[jj,] <- x_HG/sum(x_HG)
		
		csm <- csm + cj[jj]
	}
		
	##########################################################################################
	# Define the four prevalence functions
	
	if(TRUE){ # X range [-2,2]
		# mixture of normal pdf
		f1 <- function(x,mu1=-1.0,mu2=1.0,sig1=.6,sig2=.7,wt=.6){
			wt*exp(-(x-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(x-mu2)^2/(2*sig2^2))
		}
		# AEPD function
		f2 <- function(x,k=1.2,sg=2,a=2){0.7*exp((-(x>0)*k^a-(x<0)/k^a)*abs(x/sg)^a)}
		# sampled sin function
		f3 <- function(x){0.2*(sin(pi*(x-0.5)/2.5)+2)*(exp(-(x+(x+0.5)^2*(x> -0.5))/6))}
		# exponential function
		f4 <- function(x,a=0,b=2,c=4){c*exp(a+b*x)/(5+5*exp(a+b*x))}
  }
	
	if(FALSE){ # X range [-3,3]
		# mixture of normal pdf
		f1 <- function(x,mu1=-1.5,mu2=1.5,sig1=.6,sig2=.8,wt=.6){
			wt*exp(-(x-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(x-mu2)^2/(2*sig2^2))
		}
		# AEPD function
		f2 <- function(x,k=1.2,sg=2.5,a=2){0.7*exp((-(x>0)*k^a-(x<0)/k^a)*abs(x/sg)^a)}
		# sampled sin function
		f3 <- function(x){0.2*(sin(pi*(x+0.3)/2.5)+2)*(exp(-(x+(x-0.3)^2*(x>0.3))/6))}
		# exponential function
		f4 <- function(x,a=1,b=1.5,c=4){c*exp(a+b*x)/(6+6*exp(a+b*x))}
	}
	
	
	#########################################################################
	# Generate true outcome of each individual under specified models (mod=1 or 2)
	# Combine function 1 and 2
	######################################################################### model 1
	if(mod == 1){ # combine of function 1 and 2
	  
	  id_eta<- rep(NA,N)
	  int_1 <- integrate(f1,lbd_1,rbd_1)$value/(rbd_1-lbd_1)
		ff2   <- function(x,f2){f2(x)*dnorm(x)}
	  int_2 <- integrate(ff2,-15,15,f2=f2)$value
		
	  fun_1 <- f1(x_1) - int_1
	  fun_2 <- f2(x_2) - int_2
		Y_ind <- rlnorm(N,fun_1+fun_2+Wmat%*%beta,sig)
    Y_ind_HG <- Y_ind[oid]
		
	}

	########################################################################## model 2
	if(mod == 2){ # combine of function 3 and 4
	  
	  id_eta<- rep(NA,N)
	  int_1 <- integrate(f3,lbd_1,rbd_1)$value/(rbd_1-lbd_1)
		ff2   <- function(x,f2){f2(x)*dnorm(x)}
	  int_2 <- integrate(ff2,-15,15,f2=f4)$value
		
	  fun_1 <- f3(x_1) - int_1
	  fun_2 <- f4(x_2) - int_2
		Y_ind <- rlnorm(N,fun_1+fun_2+Wmat%*%beta,sig)
    Y_ind_HG <- Y_ind[oid]
	}

	#=========================================================

  #sd(Y_ind)
  #hist(Y_ind)
  #summary(lm(log(Y_ind)~Wmat[,-1]))

  # create true pooled biomarker level
  Pool_id <- rep(1:J,cj)
	Imat    <- cbind(0,Pool_id,c(t(Omat)))
	Imat_HG <- cbind(0,Pool_id,c(t(Omat_HG)))
	
	#aggregate(Imat[,3]~Imat[,2],FUN=sum )
	#Y_pool1  <- apply(matrix(Y_ind,J,cj,byrow=TRUE)*Omat,1,sum)
	
	Y_pool    <- aggregate(Y_ind*Imat[,3]~Imat[,2],FUN=sum)[,2]
	Y_pool_HG <- aggregate(Y_ind_HG*Imat_HG[,3]~Imat_HG[,2],FUN=sum)[,2]
	# sum(abs(Y_pool1-Y_pool)) # should be 0
	Zmat    <- matrix(NA,J,2+max(cj))
	Zmat_HG <- matrix(NA,J,2+max(cj))
	j1 <- 0
	for(jj in 1:J){
	  Zmat[jj,]    <- c(Y_pool[jj],cj[jj],(j1+1):(j1+cj[jj]))
		Zmat_HG[jj,] <- c(Y_pool_HG[jj],cj[jj],(j1+1):(j1+cj[jj]))
		j1 <- j1+cj[jj]
	}
  return(list(Y_ind=Y_ind,Y_pool=Y_pool,Wmat=Wmat,Omat=Omat,Imat=Imat,
	            Zmat=Zmat,x_1=x_1,x_2=x_2,
							Y_ind_HG=Y_ind_HG,Y_pool_HG=Y_pool_HG,
							Wmat_HG=Wmat_HG,Omat_HG=Omat_HG,Imat_HG=Imat_HG,
	            Zmat_HG=Zmat_HG,x_1_HG=x_1_HG,x_2_HG=x_2_HG))
}

###################################################################################
# Imat: a matrix of cj*J rows; 
#   col1=0, col2=pool id what ith individual was assigned to
#   col3=weight of ith individual in the pool 
# Zmat: a matrix of J rows
#   col1=Pooled biomarker level
#   col2=cj (pool size)
#   col3-col(2+cj)=indices of individuals assigned to jth pool
###################################################################################

# No measurement error using Gaussian Predictive Process
# cj is a vector of length J containing pool sizes
Bayes.Pool.GPP <- function(Y_pool,Wmat,Imat,Zmat,x_1,x_2,mcw_1,mcw_2,X_new_1,X_new_2,N,cj,J,
                           MCMCpara=c(1000,5,1000),beta_ini=c(1,0,0),
                           sig_ini=1,lbd,rbd,L_1,L_2){
  #MCMCpara=c(1000,5,1000);beta_ini=c(1,0,0);sig_ini=1;lbd=c(-1,-1)*2;rbd=c(1,1)*2
	nburn   <- MCMCpara[1]
  	nthin   <- MCMCpara[2]
  	nkeep   <- MCMCpara[3]
  	nsample <- nkeep*nthin
  	niter   <- nburn + nsample
	lbd_1 <- lbd[1]
	lbd_2 <- lbd[2]
	rbd_1 <- rbd[1]
	rbd_2 <- rbd[2]

	#X_new_1 <- seq(lbd_1,rbd_1,by=0.01);
  	#X_new_2 <- seq(lbd_2,rbd_2,by=0.01);
	
	KK_1 <- length(X_new_1); KK_2 <- length(X_new_2)
	K_1  <- length(x_1);   K_2  <- length(x_2)
	# ============================================================================================
	# Configuration for Gaussian Predictive Process
	#############################################################
	### determine the range of decay parameter phi for function 1
	x_length_1 = rbd_1-lbd_1
	phi_range = seq(0.01,100,by=0.01)
	lw = rep(NA,length(phi_range))
	up = lw
	for(i in 1:length(phi_range)){
	  lw[i] = matern(x_length_1*2/30,phi=phi_range[i],kappa = kap)
	  up[i] = matern(x_length_1*2/3,phi=phi_range[i],kappa = kap)
	}
	phi_min_1 = phi_range[which.min(abs(lw-0.05))]
	phi_max_1 = phi_range[which.min(abs(up-0.05))]
	#############################################################
	### determine the range of decay parameter phi for function 2
	x_length_2 = rbd_2-lbd_2
	phi_range = seq(0.01,100,by=0.01)
	lw = rep(NA,length(phi_range))
	up = lw
	for(i in 1:length(phi_range)){
	  lw[i] = matern(x_length_2*2/30,phi=phi_range[i],kappa = kap)
	  up[i] = matern(x_length_2*2/3,phi=phi_range[i],kappa = kap)
	}
	phi_min_2 = phi_range[which.min(abs(lw-0.05))]
	phi_max_2 = phi_range[which.min(abs(up-0.05))]

	logit <- function(x){log(x/(1-x))}
  phi_1_iter <- (phi_min_1+phi_max_1)/2
  trphi_1_iter <- logit((phi_1_iter-phi_min_1)/(phi_max_1-phi_min_1))
	phi_2_iter <- (phi_min_2+phi_max_2)/2
  trphi_2_iter = logit((phi_2_iter-phi_min_2)/(phi_max_2-phi_min_2))

	a_1=2;b_1=1;tau_1_iter = rgamma(1,a_1,b_1)
	a_2=2;b_2=1;tau_2_iter = rgamma(1,a_2,b_2)

	############################################################# Configure setting for function 1
	# mcw_1    = seq(lbd_1,rbd_1,length.out=L_1)   # Knots for Predictive Process 1
	knots_dist_1 = as.matrix(dist(mcw_1,mcw_1, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist_1 = matrix(NA,K_1,L_1)
	new_cross_dist_1 = matrix(NA,KK_1,L_1)
	for(j in 1:L_1){
	  cross_dist_1[,j] <- abs(x_1-mcw_1[j])
		new_cross_dist_1[,j] <- abs(X_new_1-mcw_1[j])
	}


	############################################################# Configure setting for function 2
	# mcw_2 = seq(lbd_2,rbd_2,length.out=L_2)   # Knots for Predictive Process 2
	knots_dist_2 = as.matrix(dist(mcw_2,mcw_2, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist_2 = matrix(NA,K_2,L_2)
	new_cross_dist_2 = matrix(NA,KK_2,L_2)
	for(j in 1:L_2){
	  cross_dist_2[,j] <- abs(x_2-mcw_2[j])
		new_cross_dist_2[,j] <- abs(X_new_2-mcw_2[j])
	}

	# ============================================================================================

	############################### create the matrix to store posterior samples
  fun_w_1 <- matrix(NA,L_1,nkeep); fun_w_2 <- matrix(NA,L_2,nkeep)
	fun_x_1 <- matrix(NA,K_1,nkeep); fun_x_2 <- matrix(NA,K_2,nkeep)
  fun_new_1 <- matrix(NA,KK_1,nkeep); fun_new_2 <- matrix(NA,KK_2,nkeep)
  tau_1 <- rep(NA,nkeep); tau_2 <- rep(NA,nkeep)
  phi_1 <- rep(NA,nkeep); phi_2 <- rep(NA,nkeep)
  Beta_save <- matrix(NA,p,nkeep)
	Sig_save  <- rep(NA,nkeep)
	#Er_save   <- matrix(NA,nkeep,J)
  ############################### create the distance matrix and correlation matrix for fun.1
  x.dist_row_1 = knots_dist_1[1,]
  P_w_1 = toeplitz(matern(x.dist_row_1, phi= phi_1_iter, kappa = kap)) # P star matrix
  P_w_inv_1 <- TrenchInverse(P_w_1)
  P_xw_1 = matern(cross_dist_1,phi=phi_1_iter,kappa=kap) # correlation matrix of x and x.knots
  C_star_1 = P_xw_1%*%P_w_inv_1
  ############################### create the distance matrix and correlation matrix for fun.2
  x.dist_row_2 = knots_dist_2[1,]
  P_w_2 = toeplitz(matern(x.dist_row_2, phi= phi_2_iter, kappa = kap)) # P star matrix
  P_w_inv_2 <- TrenchInverse(P_w_2)
  P_xw_2 = matern(cross_dist_2,phi=phi_2_iter,kappa=kap) # correlation matrix of x and x.knots
  C_star_2 = P_xw_2%*%P_w_inv_2

  # This code require Y_pool is a vector of pooled biomarker level
  # No measurement

  # Initialize
  #esig_iter  <- esig
  sig_iter   <- sig_ini
  beta_iter  <- beta_ini
	fun_1_iter <- rep(0,N)
	fun_2_iter <- rep(0,N)
	
	Y_iter     <- rep(Y_pool, cj)
	
  Xb_iter  <- as.numeric(Wmat%*%beta_iter)
  XW       <- as.numeric(Xb_iter + fun_1_iter + fun_2_iter)

  
  #er_iter <- rep(0.0,J)

  iter    <- 1
  ikeep   <- 1

  while(iter <= niter){

    # Sample latent individual biomarker level ======================
		if(max(cj)==1){
		  LY_iter <- log(Y_pool)
		}
		else{
      Y_iter <- Yindupdate(Y_iter,Y_pool,XW,Imat,Zmat,sig_iter,J,N)
		  LY_iter <- log( Y_iter )
		}
		
    #print(Zvec)


    #LY_iter <- log(Y_ind)

    if(iter>=0){
		
		# Update function I ======================================================
    fun_other_1_iter <- fun_2_iter
    phi_set_1 <- c(phi_1_iter,trphi_1_iter,phi_min_1,phi_max_1)
    fun.1 <- Bayes.PP.FUN(fun_other_1_iter,C_star_1,P_w_inv_1,tau_1_iter,
                              L_1,id0_mcw_1,a_1,b_1,phi_set_1,
                              x.dist_row_1,cross_dist_1,P_w_1,P_xw_1,LY_iter,Xb_iter,sig_iter)


    fun_w_1_iter <- fun.1$fun_w
    fun_1_iter <- fun.1$fun
    tau_1_iter <- fun.1$tau
    phi_1_iter <- fun.1$phi
    trphi_1_iter <- fun.1$trphi
    P_w_1 <- fun.1$P_w
    P_xw_1 <- fun.1$P_xw
    P_w_inv_1 <- fun.1$P_w_inv
    C_star_1 <- fun.1$C_star

    # Update function II =====================================================
    fun_other_2_iter <- fun_1_iter
    phi_set_2 <- c(phi_2_iter,trphi_2_iter,phi_min_2,phi_max_2)
    fun.2 <- Bayes.PP.FUN(fun_other_2_iter,C_star_2,P_w_inv_2,tau_2_iter,
                              L_2,id0_mcw_2,a_2,b_2,phi_set_2,
                              x.dist_row_2,cross_dist_2,P_w_2,P_xw_2,LY_iter,Xb_iter,sig_iter)

    fun_w_2_iter <- fun.2$fun_w
    fun_2_iter <- fun.2$fun
    tau_2_iter <- fun.2$tau
    phi_2_iter <- fun.2$phi
    trphi_2_iter <- fun.2$trphi
    P_w_2 <- fun.2$P_w
    P_xw_2 <- fun.2$P_xw
    P_w_inv_2 <- fun.2$P_w_inv
    C_star_2 <- fun.2$C_star

    }
		#fun_1_iter <- fun1
		#fun_2_iter <- fun2
		
		eta_iter <- as.numeric(LY_iter-fun_1_iter-fun_2_iter)

    # Sample regression coefficients ================================
    Cor_iter  <- solve(t(Wmat)%*%Wmat)
    BSig_iter <- Cor_iter*(sig_iter**2)
    Mu_iter   <- Cor_iter%*%colSums(Wmat*eta_iter)
		#beta_iter  <- as.numeric( mvrnorm(n = 1, mu=Mu_iter, Sigma=BSig_iter) )
    beta_iter <- (t(chol(BSig_iter)))%*%rnorm(p) + Mu_iter
    Xb_iter   <- as.numeric(Wmat%*%beta_iter)
    XW        <- as.numeric(Xb_iter + fun_1_iter + fun_2_iter)

    # Estiamte sig_iter

    sig_iter <- sqrt(1/rgamma(1,N/2+0.5,sum((eta_iter-Xb_iter)^2)/2+0.5))

    # Save and Proceed

    if((iter>nburn)&(iter%%nthin==0)){
      Beta_save[,ikeep] <- beta_iter
      Sig_save[ikeep]   <- sig_iter
			#################### save for fun.1
      fun_w_1[,ikeep] <- fun_w_1_iter
      #fun_x_1[,ikeep] <- fun_x_1_iter
      tau_1[ikeep] <- tau_1_iter
      phi_1[ikeep] <- phi_1_iter
      #################### save for fun.2
      fun_w_2[,ikeep] <- fun_w_2_iter
      #fun_x_2[,ikeep] <- fun_x_2_iter
      tau_2[ikeep] <- tau_2_iter
      phi_2[ikeep] <- phi_2_iter
      ############################################# Predict for fun.1
      PP_cross_new_1 = matern(new_cross_dist_1,phi=phi_1_iter,kappa=kap)
      CC_new_1 = PP_cross_new_1%*%P_w_inv_1
      AA_1 <- rep(1,KK_1)-diag(tcrossprod(CC_new_1,PP_cross_new_1))
      AA_1[AA_1<=0] <- 0
      fun_new_1[,ikeep] <- rnorm(KK_1,CC_new_1%*%fun_w_1_iter,AA_1/tau_1_iter)
      ############################################# Predict for fun.2
      PP_cross_new_2 = matern(new_cross_dist_2,phi=phi_2_iter,kappa=kap)
      CC_new_2 = PP_cross_new_2%*%P_w_inv_2
      AA_2 <- rep(1,KK_2)-diag(tcrossprod(CC_new_2,PP_cross_new_2))
      AA_2[AA_2<=0] <- 0
      fun_new_2[,ikeep] <- rnorm(KK_2,CC_new_2%*%fun_w_2_iter,AA_2/tau_2_iter)

			##################
			#print(ikeep)
      ikeep <- ikeep + 1
    }
    #print(iter)
    iter <- iter + 1

  }
	
  return(list(beta=Beta_save,sig=Sig_save,X_new_1=X_new_1,X_new_2=X_new_2,
	            fun_new_1=fun_new_1,fun_new_2=fun_new_2,X_1=mcw_1,X_2=mcw_2,
							tau_1=tau_1,tau_2=tau_2,phi_1=phi_1,phi_2=phi_2,
							fun_w_1=fun_w_1,fun_w_2=fun_w_2))

}


