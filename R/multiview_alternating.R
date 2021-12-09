#' Multi-View : Alternating Diffusion
#' 
#' for reusability, use markov transition kernels only. 
#' 
#' 
#' @export
mv_alternating <- function(markov1, markov2, ndim=2){
  # inputs
  if (nrow(markov1)!=nrow(markov2)){
    stop("* mv_alternating : two transition kernels have different numbers of rows.")
  }
  
  # stepsize
  markov12 = markov1%*%markov2 
  step12   = aux_entropyrule_markov(markov12)
  
  # matrix multiplication
  Pout = markov12
  for (i in 1:(step12-1)){
    Pout = markov12%*%Pout
  }
  
  
  # STEP 4. Embedding
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize   = step12
  return(output)
}