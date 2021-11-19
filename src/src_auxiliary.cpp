#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
 * (01) cpp_distance      : alternative to Rfast::Dist
 * (02) cpp_effective     : partial routine for effective resistance computation
 * (03) cpp_effective_sym : to be run inside "aux_effectivesym"
 */

// (01) cpp_distance : alternative to Rfast::Dist ==============================
// [[Rcpp::export]]
arma::mat cpp_distance(arma::mat &X){
  int N = X.n_rows;
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = arma::norm(X.row(i)-X.row(j),2);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// (02) cpp_effective : partial routine for effective resistance computation ===
// [[Rcpp::export]]
arma::mat cpp_effective(arma::mat &X){
  // parameter and output
  int N = X.n_rows;
  arma::mat output(N,N,fill::zeros);
  
  // iteration
  double theval = 0.0;
  for (int k=0;k<N;k++){
    for (int j=(k+1);j<N;j++){
      theval = X(k,k)+X(j,j)-2.0*X(k,j);
      output(k,j) = theval;
      output(j,k) = theval;
    }
  }
  
  // return
  return(output);
}

// (03) cpp_effective_sym ======================================================
// [[Rcpp::export]]
arma::mat cpp_effective_sym(arma::mat &A){
  int n = A.n_rows;
  arma::mat L    = arma::diagmat(arma::sum(A, 1))-A;
  arma::mat Linv = arma::pinv(L);
  
  arma::vec tmpvec(n,fill::zeros);
  arma::mat output(n,n,fill::zeros);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      tmpvec(i) = 1.0;
      tmpvec(j) = -1.0;
      
      output(i,j) = arma::dot(Linv*tmpvec, tmpvec);
      output(j,i) = output(i,j);
      
      tmpvec(i) = 0.0;
      tmpvec(j) = 0.0;
    }
  }
  return(output);
}
