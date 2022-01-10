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

  # construct : clustering-based transition matrix
  tmp_embed = as.matrix(phate_original(data, ndim=ndim, nbdk=nbdk, alpha=alpha, alg="mmds", potential="log")$embedding)
  aff_lab   = T4cluster::sc05Z(tmp_embed, k=round(nclust), nnbd=round(nbdk))$cluster
  
  ulabel  = unique(aff_lab)
  nulabel = length(ulabel)
  state_transition = array(0,c(nulabel, nulabel))
  for (i in 1:(N-1)){
    state_i = round(aff_lab[i])
    state_j = round(aff_lab[i+1])
    state_transition[state_i, state_j] = state_transition[state_i, state_j] + 1
  }
  
  
  
  # construct : affinity for temporal
  aff_temporal = array(0,c(N,N))
  for (i in 1:N){
    for (j in 1:N){
      state_i = round(aff_lab[i])
      state_j = round(aff_lab[j])
      aff_temporal[i,j] = state_transition[state_i, state_j]
    }
  }

  # construct : two markov matrices
  markov_data  = array(0,c(N,N))
  markov_clust = array(0,c(N,N))
  for (i in 1:N){
    tgt1 = as.vector(aff_data[i,])
    tgt2 = as.vector(aff_temporal[i,])
    
    markov_data[i,] = tgt1/base::sum(tgt1)
    markov_clust[i,]= tgt2/base::sum(tgt2)
  }
  markov_all = markov_data%*%markov_clust
  
  # optimal transition steps
  markov_step = aux_entropyrule_markov(markov_all)

  # matrix multiplication
  Pout = markov_all
  for (i in 1:(markov_step-1)){
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
  output$stepsize   = markov_step
  return(output)
}