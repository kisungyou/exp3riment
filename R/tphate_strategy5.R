#' Transition PHATE : Strategy 5
#' 
#' set 1 for same-cluster affinities
#' otherwise, use log-like decay
#' 
#' @export
tphate_strategy5 <- function(data, ndim=2, nbdk=5, alpha=2.0, scaler=0.5){
  # Inputs
  N = base::nrow(data)
  DIST_data = as.matrix(aux_dist(data))
  myc = max(0.0000001, as.double(scaler))
  
  # construct : affinity for data
  aff_tmp = aux_kernel_standard(DIST_data, round(nbdk), as.double(alpha)) 
  aff_obj = igraph::graph_from_adjacency_matrix(aff_tmp, mode="undirected", weighted=TRUE, diag=FALSE)
  aff_grp = igraph::cluster_louvain(aff_obj)
  aff_lab = as.integer(as.factor(igraph::membership(aff_grp)))
  ulabel  = unique(aff_lab)
  nulabel = length(ulabel)
  print(paste0("the number of unique states : ",nulabel))
  
  # construct : affinity for temporal
  aff_temporal = array(1,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      aff_temporal[i,j] <- aff_temporal[j,i] <- exp((log(myc)/(N-1))*(abs(i-j)))
    }
  }
  for (i in 1:nulabel){
    tgtlab = which(aff_lab==ulabel[i])
    aff_temporal[tgtlab,tgtlab]=1
  }
  
  # affinity aggregation
  affinity_machine = aff_temporal*aff_tmp
  
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
  output$membership = aff_lab
  return(output)
}
