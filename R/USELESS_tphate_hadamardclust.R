#' T-PHATE : Hadamard Product + Clustering
#' 
#' always take an (TxP) input
#' 
#' @export
tphate_hadamardclust <- function(data, ndim=2, nbdk=5, alpha=2.0, temporal=5, decay=0){
  # Inputs
  N = base::nrow(data)
  DIST_data = as.matrix(aux_dist(data))
  mydecay   = as.double(decay)
  
  # construct : affinity for data
  aff_tmp = aux_kernel_standard(DIST_data, round(nbdk), as.double(alpha)) 
  aff_obj = igraph::graph_from_adjacency_matrix(aff_tmp, mode="undirected", weighted=TRUE, diag=FALSE)
  aff_grp = igraph::cluster_louvain(aff_obj)
  aff_lab = igraph::membership(aff_grp)
  ulabel        = unique(aff_lab)
  
  affinity_data = aff_tmp
  for (i in 1:length(ulabel)){
    tgtid = which(aff_lab==ulabel[i])
    affinity_data[tgtid,tgtid] = 0
  }
  diag(affinity_data)=1

  # construct : affinity for time
  affinity_time = diag(N)
  seq1N = 1:N
  for (n in 1:N){
    now_ids = which(abs(n-seq1N)<=temporal)
    now_vec = abs(n-seq1N)[now_ids]
    affinity_time[n,now_ids] = exp(-(now_vec^mydecay))
    affinity_time[now_ids,n] = exp(-(now_vec^mydecay))
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