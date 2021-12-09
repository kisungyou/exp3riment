#' Construct PHATE-style transition kernel P=D^{-1}A
#' 
#' take (n x p) matrix or 'dist' object
#' 
#' @export
routine_transition <- function(data, nbdk=5, alpha=2.0){
  # inputs
  if (inherits(data, "dist")){
    D = data
    N = round((sqrt(8*length(D)+1)+1)/2)
  } else if (is.matrix(data)){
    N = base::nrow(data)
    D = aux_dist(data)
  }
  mat_affinity = aux_kernel_standard(D, round(nbdk), as.double(alpha))
  mat_affinity[is.na(mat_affinity)] = 1 # weird ER case : exp(-0/0tmp)
  
  output = array(0,c(N,N))
  for (i in 1:N){
    tgt = as.vector(mat_affinity[i,])
    output[i,] = tgt/base::sum(tgt)
  }
  return(output)
}