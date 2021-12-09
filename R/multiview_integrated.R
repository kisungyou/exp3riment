#' Multi-View : Integrated Diffusion
#' 
#' for reusability, use markov transition kernels only. 
#' 
#' 
#' @export
mv_integrated <- function(markov1, markov2, ndim=2){
  # inputs
  if (nrow(markov1)!=nrow(markov2)){
    stop("* mv_integrated : two transition kernels have different numbers of rows.")
  }
  
  # stepsize
  step1 = aux_entropyrule_markov(markov1)
  step2 = aux_entropyrule_markov(markov2)
  
  # matrix multiplication
  P1 = markov1
  P2 = markov2
  for (i in 1:(step1-1)){    P1 = markov1%*%P1  }
  for (j in 1:(step2-1)){    P2 = markov2%*%P2  }
  
  # embedding
  Pout = P1%*%P2
  
  # STEP 4. Embedding
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize1 = step1
  output$stepsize2 = step2
  return(output)
}