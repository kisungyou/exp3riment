#' Transition PHATE : Strategy 3
#' 
#' finite machine transition
#' 
#' @export
tphate_strategy3 <- function(data, ndim=2, nbdk=5, alpha=2.0, symmetrize=TRUE){
  # Inputs
  N = base::nrow(data)
  DIST_data = as.matrix(aux_dist(data))

  # construct : affinity for data
  aff_tmp = aux_kernel_standard(DIST_data, round(nbdk), as.double(alpha)) 
  aff_obj = igraph::graph_from_adjacency_matrix(aff_tmp, mode="undirected", weighted=TRUE, diag=FALSE)
  aff_grp = igraph::cluster_louvain(aff_obj)
  aff_lab = as.integer(as.factor(igraph::membership(aff_grp)))
  nulabel = length(unique(aff_lab))
  print(paste0("the number of unique states : ",nulabel))
  
  # construct : transition matrix
  state_transition = array(0,c(nulabel, nulabel))
  for (i in 1:(N-1)){
    state_i = round(aff_lab[i])
    state_j = round(aff_lab[i+1])
    state_transition[state_i,state_j] = state_transition[state_i,state_j] + 1
  }
  for (i in 1:nulabel){
    tgt = state_transition[i,]
    state_transition[i,] = tgt/base::sum(tgt)
  }
  
  # construct : finite-state machine
  affinity_machine = array(0,c(N,N))
  for (i in 1:N){
    for (j in 1:N){
      if (i==j){
        state_i = round(aff_lab[i])
        affinity_machine[i,i] = state_transition[state_i,state_i]
      } else {
        state_i = round(aff_lab[i])
        state_j = round(aff_lab[j])
        affinity_machine[i,j] = state_transition[state_i,state_j]
      }
    }
  }
  if (symmetrize){
    affinity_machine = (affinity_machine + t(affinity_machine))/2
  }
  affrowsum = rowSums(affinity_machine)
  affinity_machine = diag(1/sqrt(affrowsum))%*%affinity_machine%*%diag(1/sqrt(affrowsum))
  
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
  output$finite     = state_transition
  return(output)
}
