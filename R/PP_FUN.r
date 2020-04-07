
#######################################################################
# R function: This function is used to nonparametrically estimate a single unknown function of continuous
#             covariate using Full Gaussian Process.
# This function is used within BNR_GP
#
#
#
#
#
#

Bayes.PP.FUN<-function(fun_other_iter,C_star,P_w_inv,tau_iter,L,id0_mcw,a,b,phi_set,
                       x.dist_row,cross_dist,P_w,P_xw,LY_iter,Xb_iter,sig_iter){

########################################################################### update fun for fun
# fun_oter_iter is the sum of all other eat vectors corresponding to other functions
#fun_other_iter = fun_other_1_iter
#C_star   <- C_star_1
#P_w_inv  <- P_w_inv_1
#tau_iter <- tau_1_iter
#L = L_1


LXE_iter <- LY_iter-Xb_iter-fun_other_iter


# Precision matrix without sig**2
DP <- forceSymmetric(crossprod(C_star)+P_w_inv*tau_iter*sig_iter**2)
Sig_fun_iter_decomp <- chol(DP) # decomposition of presicion matrix without sig**2
Sig_fun_iter <- chol2inv(Sig_fun_iter_decomp) # covariance matrix without sig**2
mu_fun_iter <- Sig_fun_iter%*%crossprod(C_star,LXE_iter)   # mean of fun_iter; sig**2 cancelled

# sample from presicion matrix,solve uptri system
# adjust for sig_iter of LN distribution
fun_w_iter <- as.numeric(backsolve(Sig_fun_iter_decomp/sig_iter,rnorm(L))+mu_fun_iter)

#fun_w_iter <- fun_w_iter-fun_w_iter[id0_mcw]

# Let mean(fun_iter)=0


fun_iter   <- C_star%*%fun_w_iter
shift      <- sum(fun_iter)/sum(C_star)
fun_w_iter <- fun_w_iter - shift
fun_iter   <- fun_iter - mean(fun_iter)

#mf_iter  <- mean(fun_iter)
#fun_iter <- fun_iter - mf_iter
#fun_w_iter <- fun_w_iter - mf_iter
#fun_x_iter <- fun_x_iter - mf_iter

#################################### update tau_1
bb = crossprod(fun_w_iter,P_w_inv)%*%fun_w_iter/2
a_star=a+L/2;b_star=as.numeric(b+bb)
tau_iter = rgamma(1,a_star,b_star)

########################################################################### update the covariance matrix
phi_iter <- phi_set[1]
trphi_iter <- phi_set[2]
phi_min <- phi_set[3]
phi_max <- phi_set[4]

trphi_new = rnorm(1,trphi_iter,trphi_tune)
phi_new = (exp(trphi_new)*phi_max+phi_min)/(1+exp(trphi_new))

P_w_new = toeplitz(matern(x.dist_row, phi= phi_new, kappa = kap)) # Correlation matrix of x.knots, which is also a toeplitz matrix
#print("GOOD!")
P_w_inv_new = TrenchInverse(P_w_new) # Inverse the Correlation matrix of x.knots

P_xw_new = matern(cross_dist,phi=phi_new,kappa=kap) # correlation matrix of x and x.knots
C_star_new = P_xw_new%*%P_w_inv_new

#LXE_iter <- LY_iter-Xb_iter-fun_other_iter

C_star_fun_new = as.numeric(LXE_iter-C_star_new%*%fun_w_iter)
C_star_fun = as.numeric(LXE_iter-fun_iter)

rat1= -determinant(P_w_new)$modulus/2-tau_iter*crossprod(fun_w_iter,P_w_inv_new)%*%fun_w_iter/2-sum((C_star_fun_new)^2)/(2*sig_iter**2)
rat2= -determinant(P_w)$modulus/2-tau_iter*crossprod(fun_w_iter,P_w_inv)%*%fun_w_iter/2-sum((C_star_fun)^2)/(2*sig_iter**2)

rat_2 = (phi_new-phi_min)*(phi_max-phi_new)/((phi_max-phi_iter)*(phi_iter-phi_min))
rat = min(1,as.numeric(exp(rat1-rat2)*rat_2))

ber_rat = rbinom(1,1,rat)
phi_iter = ber_rat*phi_new+(1-ber_rat)*phi_iter
trphi_iter = ber_rat*trphi_new + (1-ber_rat)*trphi_iter
#accept.phi = accept.phi + ber_rat

if(ber_rat==1){ # if phi_new was accepted, then update all matrics
  P_w <- P_w_new # P star matrix
  P_w_inv <- P_w_inv_new
  P_xw <- P_xw_new # correlation matrix of x and x.knots
  C_star <- C_star_new
}

return(list("fun_w"=fun_w_iter,"fun"=fun_iter,"tau"=tau_iter,"phi"=phi_iter,
            "trphi"=trphi_iter,"P_w"=P_w,"P_w_inv"=P_w_inv,"P_xw"=P_xw,"C_star"=C_star ))
}

