#' Graphical PHATE for using Adjacency
#' 
#' Remove Landmarking Strategy at this point.
#' 
#' 
#' @export
ver2_phate_graphA <- function(input, ndim=2, nbdk=5, alpha=2.0, time_step=NULL, add_diagonal=TRUE){
  # parameters
  opt_algorithm = "mmds"
  opt_potential = "log"
  if ((is.null(time_step))&&(length(time_step)<1)){
    markov_rule = TRUE # use entropy-based determination
    markov_step = 1
  } else {
    markov_rule = FALSE
    markov_step = max(1, round(time_step))
  }
  
  # STEP 1. Use Adjacency Matrix Directly
  mat_input = aux_binarynetwork(input)
  if (add_diagonal){
    diag(mat_input) = 1
  }
  mat_P      = base::diag(1/base::rowSums(mat_input))%*%mat_input 
  
  # STEP 2. Markov Rule : Optimal Transition Steps
  if (markov_rule){
    markov_step = aux_entropyrule_markov(mat_P)
    print(paste0("* ver2_phate_graphA : optimal time step=",markov_step))
  }
  
  
  # STEP 3. Matrix Multiplication
  Pout = mat_P
  for (i in 1:(markov_step-1)){
    Pout = mat_P%*%Pout
  }
  
  # STEP 4. Embedding
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  return(output)
}