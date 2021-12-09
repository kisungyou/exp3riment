#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat spherical_merge2(arma::mat &markov1, arma::mat &markov2, double weight){
  int N = markov1.n_rows;
  
  arma::rowvec tmpvec(N,fill::zeros);
  arma::mat output(N,N,fill::zeros);
  for (int n=0; n<N; n++){
    tmpvec = (1-weight)*arma::sqrt(markov1.row(n)) + weight*arma::sqrt(markov2.row(n));
    tmpvec /= arma::norm(tmpvec, 2);
    output.row(n) = arma::square(tmpvec);
  }
  return(output);
}
