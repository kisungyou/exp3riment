#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::cube cpp_magic(arma::mat &cellxgene, int stepsize, int nbdsize, bool apply_sqrt){
  // parameters
  int n = cellxgene.n_rows;
  int p = cellxgene.n_cols;
  
  // library-size normalization
  arma::vec libsize  = arma::sum(cellxgene, 1);
  double libsize_med = arma::median(libsize);
  arma::mat data_norm(n,p,fill::zeros);
  for (int i=0; i<n; i++){
    data_norm.row(i) = (cellxgene.row(i)/libsize(i))*libsize_med;
  }
  if (apply_sqrt){
    for (int i=0; i<n; i++){
      data_norm.row(i) = arma::sqrt(data_norm.row(i));
    }
  }
  
  // construct affinity
  arma::mat pairwise_distance(n,n,fill::zeros);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      pairwise_distance(i,j) = arma::norm(data_norm.row(i)-data_norm.row(j), 2);
      pairwise_distance(j,i) = pairwise_distance(i,j);
    }
  }
  arma::rowvec vec_sorted(n,fill::zeros);
  arma::vec knn_dist(n,fill::zeros);
  for (int i=0; i<n; i++){
    vec_sorted  = arma::sort(pairwise_distance.row(i), "ascend");
    knn_dist(i) = vec_sorted(i);
  }
  double dij2 = 0.0;
  arma::mat affinity(n,n,fill::ones);
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      dij2 = std::pow(pairwise_distance(i,j), 2.0);
      affinity(i,j) = 0.5*std::exp(-dij2/(knn_dist(i)*knn_dist(i))) + 0.5*std::exp(-dij2/(knn_dist(j)*knn_dist(j)));
      affinity(j,i) = affinity(i,j);
    }
  }
  arma::mat markov(n,n,fill::zeros);
  for (int i=0; i<n; i++){
    markov.row(i) = affinity.row(i)/arma::accu(affinity.row(i));
  }
  
  // apply and stack 
  arma::cube output(n,p,stepsize+1,fill::zeros);
  output.slice(0) = data_norm;
  for (int i=1; i<=stepsize; i++){
    output.slice(i) = markov*output.slice(i-1); 
  }
  
  // scaling
  return(output);
}