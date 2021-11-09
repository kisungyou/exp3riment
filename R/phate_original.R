#' Original Version of PHATE
#' 
#' @param data matrix N observations
#' @param ndim target dimension
#' @param nbdk neighborhood size
#' @param alpha decaying kernel parameter
#' 
#' @return a list containing\describe{
#' \item{transition}{transition matrix of optimal kernel.}
#' \item{embedding}{embedding matrix.}
#' }
#' 
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
#' distEB    = stats::as.dist(fast_dist(sub_data))
#' landmark1 = phate_original(distEB, n_landmark=100)$embedding
#' landmark2 = phate_original(distEB, n_landmark=250)$embedding
#' landmark3 = phate_original(distEB, n_landmark=500)$embedding
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(landmark1, pch=19, col=sub_label, main="Case 1")
#' plot(landmark2, pch=19, col=sub_label, main="Case 2")
#' plot(landmark3, pch=19, col=sub_label, main="Case 3")
#' par(opar)
#' }
#' 
#' 
#' @export
phate_original <- function(data, ndim=2, nbdk=5, alpha=2.0, alg=c("cmds","mmds"), potential=c("log","sqrt","none"), n_landmark=100){
  # inputs
  step1_start = Sys.time()
  if (inherits(data, "dist")){
    D = data
    N = round((sqrt(8*length(D)+1)+1)/2)
  } else if (is.matrix(data)){
    N = base::nrow(data)
    D = aux_dist(data)
  } else {
    stop("* phate_original : 'dist' or 'matrix', one of two permitted.")
  }
  M = round(n_landmark) # compression down to (M x M)
  
  opt_algorithm = match.arg(alg)         # TUNABLE PARAMETER
  opt_potential = match.arg(potential)   # TUNABLE PARAMETER 
  
  fun_cmds = utils::getFromNamespace("hidden_cmds", "maotai")
  fun_mmds = utils::getFromNamespace("hidden_mmds", "maotai")
  
  time_step1 = Sys.time()-step1_start
  print(paste0("* original : ",round(as.double(time_step1),5), "secs : step 1 complete."))
  
  # build kernel matrix & markov kernel ----------------------------------------
  step2_start = Sys.time()
  
  #mat_kernel = aux_kernel_standard(D, round(nbdk), as.double(alpha))
  #mat_P      = base::diag(1/base::rowSums(mat_kernel))%*%mat_kernel 
  mats_run   = src_standard_kernel(as.matrix(D), round(nbdk), as.double(alpha))
  mat_kernel = mats_run$mat_kernel
  mat_P      = mats_run$mat_markov
  
  time_step2 = Sys.time()-step2_start
  print(paste0("* original : ",round(as.double(time_step2),5), "secs : step 2 complete."))
  
  
  if (N > M){  # reduced version
    start_landmark1 = Sys.time()
    # clustering
    pseudo_data = mat_P%*%matrix(stats::rnorm(N*50), ncol=50)
    clustering  = as.integer(factor(stats::kmeans(pseudo_data, M)$cluster))
    
    # clustering  = cluster::pam(D, M)$clustering
    
    time_landmark1 = Sys.time()-start_landmark1
    print(paste0("* original : ",round(as.double(time_landmark1),5), "secs : landmark : (1) clustering.."))
    
    
    # compression
    start_landmark2 = Sys.time()
    compression = aux_compression(mat_P, mat_kernel, M, clustering)
    P_MN = compression$P_MN
    P_NM = compression$P_NM
    P_MM = P_MN%*%P_NM
    
    time_landmark2 = Sys.time()-start_landmark2
    print(paste0("* original : ",round(as.double(time_landmark2),5), "secs : landmark : (2) compression.."))
    
    
    # optimal t
    start_landmark3 = Sys.time()
    opt.t = aux_entropyrule_markov(P_MM)
    
    time_landmark3 = Sys.time()-start_landmark3
    print(paste0("* original : ",round(as.double(time_landmark3),5), "secs : landmark : (3) entropy decision.."))
    
    # embedding
    start_landmark4 = Sys.time()
    Pout  = P_MM
    for (i in 1:(opt.t-1)){
      Pout = P_MM%*%Pout
    }
    Ylandmark = phate_original_embedding(Pout, ndim, opt_algorithm, opt_potential)
    Y = P_NM%*%Ylandmark
    
    time_landmark4 = Sys.time()-start_landmark4
    print(paste0("* original : ",round(as.double(time_landmark4),5), "secs : landmark : (4) embedding.."))
    
    # just a quick test
    # fin_P = P_NM%*%Pout
    # fin2d = phate_original_embedding(fin_P, ndim, opt_algorithm, opt_potential)
    
    # return
    output = list()
    output$transition = Pout
    output$embedding  = fin2d #Y
    return(output)
  } else {     # naive version
    # optimal t
    opt.t = aux_entropyrule(mat_kernel)
    print(paste0("optimal t=",opt.t))
    
    # embedding
    Pout  = mat_P
    for (i in 1:(opt.t-1)){
      Pout = mat_P%*%Pout
    }
    Y = phate_original_embedding(Pout, ndim, opt_algorithm, opt_potential)
    
    # return
    output = list()
    output$transition = Pout
    output$embedding  = Y
    return(output)
  }
}


# Embedding for Original PHATE --------------------------------------------
#' @keywords internal
#' @noRd
phate_original_embedding <- function(P, n_dim, algorithm, potential){
  # potential selection and transformation
  if (all(potential=="log")){
    rowMat = base::log(P+(1e-8))
  } else if (all(potential=="sqrt")){
    rowMat = base::sqrt(P)
  } else {
    rowMat = P
  }
  
  if (all(algorithm=="cmds")){
    # Classical MDS 
    # fun_embed = utils::getFromNamespace("hidden_cmds", "maotai")
    # return(fun_embed(transDist, ndim = n_dim)$embed)
    return(aux_fastpca(rowMat, ndim=n_dim))
  } else {
    # Metric MDS
    transDist = aux_dist(rowMat)
    fun_embed = utils::getFromNamespace("hidden_mmds", "maotai")
    return(fun_embed(transDist, ndim = n_dim))
  }
}