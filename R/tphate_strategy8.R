#' Transition PHATE : Strategy 8
#' 
#' alternating diffusion : data + cluster-based transition
#' 
#' @export
tphate_strategy8 <- function(data, ndim=2, nbdk=5, alpha=2.0, nclust=10){
  # Inputs
  N = base::nrow(data)
  DIST_data = as.matrix(aux_dist(data))
  nclust    = round(nclust)
  
  # construct : affinity for data
  aff_data = aux_kernel_standard(DIST_data, round(nbdk), as.double(alpha)) 
  

  # construct : clustering-based transition
  tmp_embed = phate_original(data, ndim=ndim, nbdk=nbdk, alpha=alpha, alg="mmds", potential="log")$embedding
  tmp_label = T4cluster::sc05Z(tmp_embed)$cluster
    
  
  # construct : affinity for time
  aff_temporal = array(1,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      aff_temporal[i,j] <- aff_temporal[j,i] <- exp((log(myc)/(N-1))*(abs(i-j)))
    }
  }
  
  # affinity aggregation
  affinity_machine = aff_temporal*aff_data
  
  # construct : markov transition kernel
  markov_hadamard = array(0,c(N,N))
  for (i in 1:N){
    tgt = as.vector(affinity_machine[i,])
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