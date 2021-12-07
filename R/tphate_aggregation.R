#' T-PHATE : Distance Aggregation
#' 
#' always take an (TxP) input
#' 
#' @export
tphate_aggregation <- function(data, ndim=2, nbdk=5, alpha=2.0, temporal=5, weight=0.5){
  # Inputs
  N = base::nrow(data)
  DIST_data = aux_dist(data)
  DIST_time = array(0,c(N,N))
  
  seq1N = 1:N
  for (n in 1:N){
    now_vec = abs(n-seq1N)
    DIST_time[n,] = now_vec
    DIST_time[,n] = now_vec
  }
  
  # Aggregate Distances
  DIST_all = phate_mview_ver1_merge(stats::as.dist(DIST_data), stats::as.dist(DIST_time), weight)
  
  # construct affinity
  affinity_all = aux_kernel_standard(DIST_all, round(nbdk), as.double(alpha)) 
  
  # construct : markov transition kernel
  markov_all = array(0,c(N,N))
  for (i in 1:N){
    tgt = as.vector(affinity_all[i,])
    markov_all[i,] = tgt/base::sum(tgt)
  }
  
  # optimal transition steps
  step_aggregation = aux_entropyrule_markov(markov_all)
  
  # matrix multiplication
  Pout = markov_all
  for (i in 1:(step_aggregation-1)){
    Pout = markov_all%*%Pout
  }
  
  # STEP 4. Embedding
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize   = step_aggregation
  return(output)
}