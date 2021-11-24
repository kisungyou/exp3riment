#' Original PHATE Simplified
#' 
#' 
#' 
#' @export
ver2_phate_original <- function(data, ndim=2, nbdk=5, alpha=2.0, time_step=NULL){
  # inputs
  if (inherits(data, "dist")){
    D = data
    N = round((sqrt(8*length(D)+1)+1)/2)
  } else if (is.matrix(data)){
    N = base::nrow(data)
    D = aux_dist(data)
  }
  opt_algorithm = "mmds"
  opt_potential = "log"
  if ((is.null(time_step))&&(length(time_step)<1)){
    markov_rule = TRUE # use entropy-based determination
    markov_step = 1
  } else {
    markov_rule = FALSE
    markov_step = max(1, round(time_step))
  }
  
  # STEP 1. build kernel
  mat_kernel = aux_kernel_standard(D, round(nbdk), as.double(alpha))
  mat_kernel[is.na(mat_kernel)] = 1 # weird ER case : exp(-0/0tmp)
  mat_P      = base::diag(1/base::rowSums(mat_kernel))%*%mat_kernel 
  
  # STEP 2. Markov Rule : Optimal Transition Steps
  if (markov_rule){
    markov_step = aux_entropyrule_markov(mat_P)
    print(paste0("* ver2_phate_original : optimal time step=",markov_step))
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