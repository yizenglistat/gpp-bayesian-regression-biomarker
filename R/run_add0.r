rm(list=ls(all=TRUE))

Cluster=TRUE

if(Cluster==TRUE)
{
  .libPaths(new="/home/wang528/R")
  setwd("/home/wang528/PoolAPLM")
  
  options(echo=TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  vid <- as.numeric(args)
}

if(Cluster==FALSE)
{
  vid<-1
}

library(MASS);
library(Matrix);
library(mvtnorm);
library(truncnorm);
library(coda);
library(geoR);
library(ltsa);
library(Rcpp);
library(RcppTN);

source("Support_functions_add_weighted_general.R")
source("PP_FUN.R")
Rcpp::sourceCpp("Yindupdate.cpp")



weight<- TRUE
J     <- 500            # number of pools

# ===============================================
nsim  <- 12            # number of simulations per core
# ===============================================

# pj  # pool size; 1,2,4,6,8
# mod # model selection: 1 or 2
for(copula.rho0 in 0)
{
  if(copula.rho0==0){copula.rho=0}
  if(copula.rho0==1){copula.rho=0.8}
  if(copula.rho0==2){copula.rho=-0.8}
  if(copula.rho0==3){copula.rho=0.5}
  if(copula.rho0==4){copula.rho=-0.5}
  if(copula.rho0==5){copula.rho=0.2}
  if(copula.rho0==6){copula.rho=-0.2}
  
  for(pj in c(1,8,6,4,2))
  {
    for(mod in 1:2)
    {
    cj    <- rep(pj,J)      # vector of pool size
    N     <- sum(cj)        # total number of subjects
    L     <- 2              # frequency of repeated measurements
    esig  <- 0.1            # sd of measurement error
    sig   <- 0.5            # sd of lognormal distribution
    beta  <- c(0.5,.5,-.5)
    p     <- length(beta)
    kap   <- 2.0            # smooth parameter in matern function
    Q     <- 2              # number of unknown functions
    L_1   <- 50             # number of knots for Gaussian predictive process associated with function 1 (GPP_1)
    L_2   <- 50             # number of knots for Gaussian predictive process associated with function 2 (GPP_2)
    lbd   <- rep(-1,2)*2
    rbd   <- rep(1,2)*2
    lbd_1 <- lbd[1];rbd_1 <- rbd[1]  # left and right endpoints for support of function 1
    lbd_2 <- lbd[2];rbd_2 <- rbd[2]  # left and right endpoints for support of function 2
    X_new_1 <- seq(lbd_1,rbd_1,by=0.02);
    X_new_2 <- seq(lbd_2,rbd_2,by=0.02);
    
    trphi_tune <- 0.1       # tune parameter for decay parameter phi

    beta_mat <- matrix(NA,nsim,2*p)
    beta_lb  <- matrix(NA,nsim,2*p)
    beta_ub  <- matrix(NA,nsim,2*p)
    beta_sd  <- matrix(NA,nsim,2*p)
    sig_vec  <- matrix(NA,nsim,2)
    sig_lb   <- matrix(NA,nsim,2)
    sig_ub   <- matrix(NA,nsim,2)
    sig_sd   <- matrix(NA,nsim,2)
    
    nter <- 1
    
    while(nter <= nsim){
      #########################################################################
      
      dat <- data_sim(weight=weight,mod=mod,beta,sig,esig,N,cj,J,lbd=lbd,rbd=rbd,L_1,L_2,copula.rho)
      
      # ===================================
      Y_ind  = dat$Y_ind
      Y_pool = dat$Y_pool
      Wmat   = dat$Wmat
      Omat   = dat$Omat
      Imat   = dat$Imat
      Zmat   = dat$Zmat
      x_1  = dat$x_1
      x_2  = dat$x_2
      
      Y_ind_HG  = dat$Y_ind_HG
      Y_pool_HG = dat$Y_pool_HG
      Wmat_HG   = dat$Wmat_HG
      Omat_HG   = dat$Omat_HG
      Imat_HG   = dat$Imat_HG
      Zmat_HG   = dat$Zmat_HG
      x_1_HG  = dat$x_1_HG
      x_2_HG  = dat$x_2_HG
      
      # ===================================
      
      mcw_1   <- seq(lbd_1,rbd_1,length.out=L_1)                             # Knots for Predictive Process 1
      lbd_2   <- quantile(x_2,0.01)
      rbd_2   <- quantile(x_2,0.99)
      lbd[2]  <- lbd_2
      rbd[2]  <- rbd_2
      mcw_2   <- seq(lbd_2,rbd_2,length.out=L_1)   # Knots for Predictive Process 2
      
      t0 = Sys.time()
      
      fitGPP <- Bayes.Pool.GPP(Y_pool,Wmat,Imat,Zmat,x_1,x_2,mcw_1,mcw_2,
                               X_new_1,X_new_2,N,cj,J,MCMCpara=c(5000,10,2000),beta_ini=c(1,0,0),
                               sig_ini=1,lbd=lbd,rbd=rbd,L_1=L_1,L_2=L_2)
      t1 = Sys.time()
      if(pj==1){
			  resGPP <- fitGPP
			}else{
			resGPP <- Bayes.Pool.GPP(Y_pool_HG,Wmat_HG,Imat_HG,Zmat_HG,x_1_HG,x_2_HG,mcw_1,mcw_2,
                               X_new_1,X_new_2,N,cj,J,MCMCpara=c(5000,10,2000),beta_ini=c(1,0,0),
                               sig_ini=1,lbd=lbd,rbd=rbd,L_1=L_1,L_2=L_2)
			}
      
      t2 = Sys.time()
      
      tm1 <- difftime(t1,t0,units="secs")
      tm2 <- difftime(t2,t1,units="secs")
      
      if(FALSE){
        HPDinterval( as.mcmc(t(fitGPP$beta)) )
        HPDinterval( as.mcmc(t(resGPP$beta)) )
        apply(fitGPP$beta,1,median)
        quantile(fitGPP$beta[3,],prob=c(0.025,0.975))
        
        windows()
        par(mfrow=c(2,2))
        plot(fitGPP$beta[1,],type='l')
        lines(resGPP$beta[1,],lty=2,col="red")
        plot(fitGPP$beta[2,],type='l')
        lines(resGPP$beta[2,],lty=2,col="red")
        plot(fitGPP$beta[3,],type='l')
        lines(resGPP$beta[3,],lty=2,col="red")
        
        HPDinterval( as.mcmc(t(fitGPP$beta)) )
        quantile(fitGPP$beta[3,],prob=c(0.025,0.975))
        # AEPD function
        f2 <- function(x,k=1.2,sg=2,a=2){0.7*exp((-(x>0)*k^a-(x<0)/k^a)*abs(x/sg)^a)}
        # mixture of normal pdf
        f1 <- function(x,mu1=-1.0,mu2=1.0,sig1=.6,sig2=.7,wt=.6){
          wt*exp(-(x-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(x-mu2)^2/(2*sig2^2))
        }
        # sampled sin function
        f3 <- function(x){0.2*(sin(pi*(x-0.5)/2.5)+2)*(exp(-(x+(x+0.5)^2*(x> -0.5))/6))}
        # exponential function
        f4 <- function(x,a=0,b=2,c=4){c*exp(a+b*x)/(5+5*exp(a+b*x))}
        int_1 <- integrate(f1,lbd_1,rbd_1)$value/(rbd_1-lbd_1)
        ff2   <- function(x,f2){f2(x)*dnorm(x)}
        int_2 <- integrate(ff2,lbd_2,rbd_2,f2=f2)$value
        
        par(mfrow=c(1,2))
        plot(X_new_1,f1(X_new_1)-int_1,type='l',ylim=c(-0.7,0.5))
        lines(X_new_1,apply(fitGPP$fun_new_1,1,mean),col='blue')
        lines(X_new_1,apply(resGPP$fun_new_1,1,mean),col='red')
        
        plot(X_new_2,f2(X_new_2)-int_2,type='l',ylim=c(-0.7,0.5))
        lines(X_new_2,apply(fitGPP$fun_new_2,1,mean),col='blue')
        lines(X_new_2,apply(resGPP$fun_new_2,1,mean),col='red')
        
      }
      
      bqt <- apply(rbind(fitGPP$beta,resGPP$beta),1,quantile,prob=c(0.5,0.025,0.975))
      bsd <- apply(rbind(fitGPP$beta,resGPP$beta),1,sd)
      sqt <- apply(rbind(fitGPP$sig, resGPP$sig) ,1,quantile,prob=c(0.5,0.025,0.975))
      ssd <- apply(rbind(fitGPP$sig, resGPP$sig) ,1,sd)
      
      RME_FG1 <- apply(fitGPP$fun_new_1,1,mean)
      RME_FG2 <- apply(fitGPP$fun_new_2,1,mean)
      HME_FG1 <- apply(resGPP$fun_new_1,1,mean)
      HME_FG2 <- apply(resGPP$fun_new_2,1,mean)
      
      RMD_FG1 <- apply(fitGPP$fun_new_1,1,median)
      RMD_FG2 <- apply(fitGPP$fun_new_2,1,median)
      HMD_FG1 <- apply(resGPP$fun_new_1,1,median)
      HMD_FG2 <- apply(resGPP$fun_new_2,1,median)
      
      
      ff1  <- c(X_new_1,RME_FG1,RMD_FG1,HME_FG1,HMD_FG1)
      ff2  <- c(X_new_2,RME_FG2,RMD_FG2,HME_FG2,HMD_FG2)
      
      # Save results to .txt file =========================================================================
      write(c(t(bqt),bsd,t(sqt),ssd),file=paste0("Mod",mod,"_cop",copula.rho0,"_J",J,"_cj",pj,"_D",vid,"_",weight,"_Par.txt"),append=TRUE)
      write(ff1,file=paste0("Mod",mod,"_cop",copula.rho0,"_J",J,"_cj",pj,"_D",vid,"_",weight,"_Fun1.txt"),append=TRUE)
      write(ff2,file=paste0("Mod",mod,"_cop",copula.rho0,"_J",J,"_cj",pj,"_D",vid,"_",weight,"_Fun2.txt"),append=TRUE)
      write(c(tm1,tm2),file=paste0("Mod",mod,"_cop",copula.rho0,"_J",J,"_cj",pj,"_D",vid,"_",weight,"_Time.txt"),append=TRUE)
      
      
      print(paste("Finish simulation ------ ",nter))
      nter <- nter + 1
      }
    }
  }
}


