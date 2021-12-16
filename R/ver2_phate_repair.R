#' Original PHATE Simplified + Metric Repair
#' 
#' @examples 
#' \dontrun{
#' # load the data
#' data(Embryoid, package="exp3riment")
#' 
#' nsample = 2000
#' sub_id  = sample(1:nrow(Embryoid$data), nsample)
#' 
#' sub_data  = Embryoid$data[sub_id,]
#' sub_label = Embryoid$label[sub_id]
#' 
#' landmark1 = ver2_phate_repair(sub_data)$embedding
#' landmark2 = ver2_phate_repair(sub_data)$embedding
#' landmark3 = ver2_phate_repair(sub_data)$embedding
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(landmark1, pch=19, col=sub_label, main="Case 1")
#' plot(landmark2, pch=19, col=sub_label, main="Case 2")
#' plot(landmark3, pch=19, col=sub_label, main="Case 3")
#' par(opar)
#' }
#' 
#' @export
ver2_phate_repair <- function(data, ndim=2, nbdk=5, alpha=2.0, time_step=NULL){
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
  Y = phate_repair_embedding(Pout, round(ndim), opt_algorithm, opt_potential)
  
  # Return
  output = list()
  output$transition = Pout
  output$embedding  = Y
  return(output)
}


# apply repair embedding --------------------------------------------------
#' @keywords internal
#' @noRd
phate_repair_embedding <- function(P, n_dim, algorithm, potential){
  # potential selection and transformation
  if (all(potential=="log")){
    rowMat = base::log(P+(1e-8))
  } else if (all(potential=="sqrt")){
    rowMat = base::sqrt(P)
  } else {
    rowMat = P
  }
  
  dist_mat    = as.matrix(aux_dist(rowMat))
  dist_repair = stats::as.dist(increase_only_metric_repair(dist_mat))
  
  if (all(algorithm=="cmds")){
    # Classical MDS 
    fun_embed = utils::getFromNamespace("hidden_cmds", "maotai")
    return(fun_embed(dist_repair, ndim = n_dim)$embed)
  } else {
    # Metric MDS
    fun_embed = utils::getFromNamespace("hidden_mmds", "maotai")
    return(fun_embed(dist_repair, ndim = n_dim))
  }
}