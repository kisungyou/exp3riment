#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List src_standard_kernel(arma::mat& D, arma::uword nbdk, double alpha){
  arma::uword N = D.n_rows;
  
  // nearest distance
  arma::rowvec rowvecD(N,fill::zeros);
  arma::vec near_distance(N,fill::zeros);
  for (arma::uword n=0; n<N; n++){
    rowvecD = arma::sort(D.row(n), "ascend");
    near_distance(n) = rowvecD(nbdk);
  }
  
  // construct a kernel
  double term1 = 0.0;
  double term2 = 0.0;
  arma::mat output_kernel(N,N,fill::ones);
  for (arma::uword i=0; i<(N-1); i++){
    for (arma::uword j=(i+1); j<N; j++){
      term1 = std::exp(-std::pow(-(D(i,j)/near_distance(i)), alpha));
      term2 = std::exp(-std::pow(-(D(i,j)/near_distance(j)), alpha));
      output_kernel(i,j) = 0.5*(term1+term2);
      output_kernel(j,i) = output_kernel(i,j);
    }
  }
  
  // construct a row-stochastic matrix 
  arma::mat output_markov(N,N,fill::zeros);
  for (arma::uword n=0; n<N; n++){
    rowvecD = output_kernel.row(n);
    output_markov.row(n) = rowvecD/arma::accu(rowvecD);
  }
  
  // return
  return(Rcpp::List::create(Rcpp::Named("mat_kernel")=output_kernel,
                            Rcpp::Named("mat_markov")=output_markov));
}