#include <Rcpp.h>
#include <RcppTN.h>
// [[Rcpp::depends(RcppTN)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector sampleLatentBM(NumericVector Y_ind,NumericVector Y_pool,NumericVector XW,
                         NumericMatrix Imat, NumericMatrix Zmat,
                         double sig_iter,int J, int N) {
  int i,j,k,l,i1,ij,cj;
  NumericVector Y(N);
  double Yj,w1,y1,x1,wi,yi,xi,wk,yk,y1new,yinew;
  double EUB,UB,r11,r12,r21,r22,r31,r32,r4,rr,lu,accept;
  Y = Y_ind;
  for(j=0;j<J;j++){   // update by group
    Yj = Y_pool[j];   // biomarker level of jth pool
    cj = Zmat(j,1);   // size of jth pool; need to upate cj-1 items
    i1 = Zmat(j,2)-1; // id of first ind. in the jth pool; -1
    w1 = Imat(i1,2);  // weight of first ind.
    y1 = Y[i1];       // biomarker level of 1st ind.
    x1 = XW[i1];      // predictor of 1st ind.

    for(l=1;l<cj;l++){
      i   = Zmat(j,2+l)-1; // id of ind. to be updated; -1
      wi  = Imat(i,2);     // weight of ith ind.
      yi  = Y[i];
      xi  = XW[i];         // predictor of ith ind.
      EUB = Yj/w1;
      for(ij=1;ij<cj;ij++){
        if(ij!=l){
          k  = Zmat(j,2+ij)-1;
          wk = Imat(k,2);
          yk = Y[k];
          EUB -=  yk*wk/w1 ;
        }
      }
      UB    = log(EUB*w1/wi); // upper bound of TN
      yinew  = exp( RcppTN::rtn1(log(yi), sig_iter, -HUGE_VAL, UB) );
      y1new = EUB - yinew*wi/w1;
      r11 = R::dlnorm(y1new,x1,sig_iter,TRUE);
      r12 = R::dlnorm(yinew,xi,sig_iter,TRUE);
      r21 = R::dlnorm(y1,x1,sig_iter,TRUE);
      r22 = R::dlnorm(yi,xi,sig_iter,TRUE);
      r31 = R::pnorm((UB-log(yi))/sig_iter,0,1,TRUE,TRUE);
      r32 = R::pnorm((UB-log(yinew))/sig_iter,0,1,TRUE,TRUE);
      r4 = log(yinew) - log(yi);
      rr = r11+r12-r21-r22+r31-r32+r4;
      lu = log(R::runif(0,1));
      accept = 0;
      if(rr>lu){
        accept=1;
      }
      yi = yinew*accept + yi*(1-accept);
      y1 = y1new*accept + y1*(1-accept);

      Y[i1] = y1;
      Y[i]  = yi;
    }


  }

  return Y;
}



//RcppTN::rtn1(g1j/G1j, 1/sqrt(G1j), 0, HUGE_VAL)

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//



