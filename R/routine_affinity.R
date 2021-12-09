#' Construct PHATE-style affinity
#' 
#' take (n x p) matrix or 'dist' object
#' 
#' @export
routine_affinity <- function(data, nbdk=5, alpha=2.0){
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
  return(mat_affinity)
}