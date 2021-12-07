#' T-PHATE : Integrated Diffusion
#' 
#' always take an (TxP) input
#' @export
tphate_hadamard <- function(data, ndim=2, nbdk=5, alpha=2.0, temporal=5){
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
  affinity_hadamard = affinity_data*affinity_time
  
  # construct : markov transition kernel
  markov_hadamard = array(0,c(N,N))
  for (i in 1:N){
    tgt = as.vector(affinity_hadamard[i,])
    markov_hadamard[i,] = tgt/base::sum(tgt)
  }
  
  # optimal transition steps
  step_hadamard = aux_entropyrule_markov(markov_hadamard)
  
  # matrix multiplication
  Pout = markov_hadamard
  for (i in 1:(step_hadamard-1)){
    Pout = markov_hadamard%*%Pout
  }
  
  # STEP 4. Embedding
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize   = step_hadamard
  return(output)
}