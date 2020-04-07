#include <Rcpp.h>
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

// See McMahan's MHS_GTR_IP file for supplemental information
// N = number of individuals
// np = number of pools current individual is in
// p = current vector of probabilities, i.e. inverse link function
// Y = matrix with pool information for individuals
// Z = matrix with observations and other information

NumericVector createXbeta(int N,int nbeta, NumericVector Xbeta, NumericMatrix X,NumericMatrix beta) {	

	for(int i=0; i<N; i++){
		for(int j=0; j<nbeta; j++){
			Xbeta(i) += (X(i,j) * beta(i,j));
		}
	}
	return Xbeta;
}