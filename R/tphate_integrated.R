#' T-PHATE : Integrated Diffusion
#' 
#' always take an (TxP) input
#' @export
tphate_integrated <- function(data, ndim=2, nbdk=5, alpha=2.0, temporal=5){
  # Inputs
  N = base::nrow(data)
  DIST_data = aux_dist(data)
  
  # construct : affinity
  affinity_data = aux_kernel_standard(DIST_data, round(nbdk), as.double(alpha)) 
  affinity_time = diag(N)
  seq1N = 1:N
  for (n in 1:N){
    now_ids = which(abs(n-seq1N)<=temporal)
    now_vec = abs(n-seq1N)[now_ids]
    affinity_time[n,now_ids] = exp(-now_vec)
    affinity_time[now_ids,n] = exp(-now_vec)
  }
  
  # construct : markov transition kernel
  markov_data = array(0,c(N,N))
  markov_time = array(0,c(N,N))
  for (i in 1:N){
    tgt = as.vector(affinity_data[i,])
    markov_data[i,] = tgt/base::sum(tgt)
  }
  for (i in 1:N){
    tgt = as.vector(affinity_time[i,])
    markov_time[i,] = tgt/base::sum(tgt)
  }
  
  # optimal transition steps
  step_data = aux_entropyrule_markov(markov_data)
  step_time = aux_entropyrule_markov(markov_time)
  
  # matrix multiplication
  P_data = markov_data
  P_time = markov_time
  for (i in 1:(step_data-1)){
    P_data = markov_data%*%P_data
  }
  for (j in 1:(step_time-1)){
    P_time = markov_time%*%P_time
  }
  Pout = P_data%*%P_time
  
  # STEP 4. Embedding
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize_data = step_data
  output$stepsize_time = step_time
  return(output)
}