#' Multi-View : Spherically Aggregated Diffusion
#' 
#' for reusability, use markov transition kernels only. 
#' 
#' 
#' @export
mv_spherical <- function(markov1, markov2, ndim=2, weight=0.5){
  # inputs
  if (nrow(markov1)!=nrow(markov2)){
    stop("* mv_spherical : two transition kernels have different numbers of rows.")
  }
  
  # stepsize
  step1 = aux_entropyrule_markov(markov1)
  step2 = aux_entropyrule_markov(markov2)
  
  # matrix multiplication
  P1 = markov1
  P2 = markov2
  for (i in 1:(step1-1)){    P1 = markov1%*%P1  }
  for (j in 1:(step2-1)){    P2 = markov2%*%P2  }
  
  # aggregation
  if (weight <= 0){
    Pout = P1
  } else if (weight >= 1){
    Pout = P2
  } else {
    Pout = spherical_merge2(P1, P2, weight)
  }
  
  # run
  opt_algorithm = "mmds"
  opt_potential = "log"
  Y = phate_original_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  output$stepsize1 = step1
  output$stepsize2 = step2
  return(output)
}

# data("Digits")
# data1 = Digits$fou
# data2 = Digits$kar
# 
# markov1 = routine_transition(data1)
# markov2 = routine_transition(data2)
# 
# 
# mv_alt = mv_alternating(markov1, markov2)$embedding
# mv_int = mv_integrated(markov1, markov2)$embedding
# mv_sph = mv_spherical(markov1, markov2)$embedding
# 
# 
# dlab = as.factor(Digits$label)
# par(mfrow=c(1,3))
# plot(mv_alt, pch=19, col=dlab, main="Alternating")
# plot(mv_int, pch=19, col=dlab, main="Integrated")
# plot(mv_sph, pch=19, col=dlab, main="Spherical")