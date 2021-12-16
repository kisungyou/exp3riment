#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat increase_only_metric_repair(arma::mat &D){
  int N = D.n_rows;
  int counter = 0;
  
  arma::mat output = D;
  
  for (int k=0; k<N; k++){
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        if (output(i,j) > output(i,k) + output(k,j)){
          counter += 1;
          output(i,j) = output(i,k) + output(k,j);
          output(j,i) = output(i,j);
        }
      }
    }
  }
  Rcpp::Rcout << "* metric_repair : violation numbers=" << counter << std::endl;
  return(output);
}