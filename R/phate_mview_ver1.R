#' MultiView PHATE Version 1 : Distance Merging (2-View Only)
#' 
#' weight=0 means using dist1 only.
#' 
#' @return a list containing\describe{
#' \item{transition}{transition matrix of optimal kernel.}
#' \item{embedding}{embedding matrix.}
#' }
#' 
#' @examples 
#' \dontrun{
#' data(iris)
#' A = as.matrix(iris[,1:4]) + 0.1*matrix(rnorm(150*4), ncol=4)
#' B = as.matrix(iris[,1:4]) + 0.1*matrix(rnorm(150*4), ncol=4)
#' 
#' D1 = stats::dist(A)
#' D2 = stats::dist(B)
#' output = phate_mview_ver1(D1, D2, weight=0.5)
#' }
#' 
#' 
#' 
#' 
#' @export
phate_mview_ver1 <- function(dist1, dist2, weight=0.5, nbdk=5, ndim=2, alpha=2.0, alg=c("cmds","mmds"), potential=c("log","sqrt","none"), n_landmark=100){
  # inputs
  step1_start = Sys.time()
  if (!inherits(dist1, "dist")){
    stop("* phate_mview_ver1 : 'dist1' is not 'dist' object.")
  }
  if (!inherits(dist2, "dist")){
    stop("* phate_mview_ver1 : 'dist2' is not 'dist' object.")
  }
  N1 = round((sqrt(8*length(dist1)+1)+1)/2)
  N2 = round((sqrt(8*length(dist2)+1)+1)/2)
  if (N1!=N2){
    stop("* phate_mview_ver1 : two input distances are not matching.")
  }
  weight = as.double(weight)
  if ((weight < 0)||(weight > 1)){
    stop("* phate_mview_ver1 : weight in [0,1].")
  }
  D = phate_mview_ver1_merge(dist1, dist2, weight)
  N = base::nrow(D)
  M = round(n_landmark) # compression down to (M x M)
  
  opt_algorithm = match.arg(alg)         # TUNABLE PARAMETER
  opt_potential = match.arg(potential)   # TUNABLE PARAMETER 
  
  fun_cmds = utils::getFromNamespace("hidden_cmds", "maotai")
  fun_mmds = utils::getFromNamespace("hidden_mmds", "maotai")
  
  time_step1 = Sys.time()-step1_start
  print(paste0("* mview_ver1 : ",round(as.double(time_step1),5), "secs : step 1 complete."))
  
  # build kernel matrix & markov kernel ----------------------------------------
  step2_start = Sys.time()
  
  #mat_kernel = aux_kernel_standard(D, round(nbdk), as.double(alpha))
  #mat_P      = base::diag(1/base::rowSums(mat_kernel))%*%mat_kernel 
  mats_run   = src_standard_kernel(as.matrix(D), round(nbdk), as.double(alpha))
  mat_kernel = mats_run$mat_kernel
  mat_P      = mats_run$mat_markov
  
  time_step2 = Sys.time()-step2_start
  print(paste0("* mview_ver1 : ",round(as.double(time_step2),5), "secs : step 2 complete."))
  
  
  if (N > M){  # reduced version
    start_landmark1 = Sys.time()
    # clustering
    pseudo_data = mat_P%*%matrix(stats::rnorm(N*50), ncol=50)
    clustering  = as.integer(factor(stats::kmeans(pseudo_data, M)$cluster))
    
    # clustering  = cluster::pam(D, M)$clustering
    
    time_landmark1 = Sys.time()-start_landmark1
    print(paste0("* mview_ver1 : ",round(as.double(time_landmark1),5), "secs : landmark : (1) clustering.."))
    
    
    # compression
    start_landmark2 = Sys.time()
    compression = aux_compression(mat_P, mat_kernel, M, clustering)
    P_MN = compression$P_MN
    P_NM = compression$P_NM
    P_MM = P_MN%*%P_NM
    
    time_landmark2 = Sys.time()-start_landmark2
    print(paste0("* mview_ver1 : ",round(as.double(time_landmark2),5), "secs : landmark : (2) compression.."))
    
    
    # optimal t
    start_landmark3 = Sys.time()
    opt.t = aux_entropyrule_markov(P_MM)
    
    time_landmark3 = Sys.time()-start_landmark3
    print(paste0("* mview_ver1 : ",round(as.double(time_landmark3),5), "secs : landmark : (3) entropy decision.."))
    
    # embedding
    start_landmark4 = Sys.time()
    Pout  = P_MM
    for (i in 1:(opt.t-1)){
      Pout = P_MM%*%Pout
    }
    Ylandmark = phate_original_embedding(Pout, ndim, opt_algorithm, opt_potential)
    Y = P_NM%*%Ylandmark
    
    time_landmark4 = Sys.time()-start_landmark4
    print(paste0("* mview_ver1 : ",round(as.double(time_landmark4),5), "secs : landmark : (4) embedding.."))
    
    # return
    output = list()
    output$transition = Pout
    output$embedding  = Y
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



# merge two distances -----------------------------------------------------
#' @keywords internal
#' @noRd
phate_mview_ver1_merge <- function(dist1, dist2, weight){
  D1 = as.matrix(dist1)/base::mean(dist1) # normalization 
  D2 = as.matrix(dist2)/base::mean(dist2) # normalization
  N  = base::nrow(D1)
  t  = as.double(weight)
  
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      output[i,j] <- sqrt(((1-t)*((D1[i,j])^2)) + (t*((D2[i,j])^2)))
      output[j,i] <- output[i,j]
    }
  }
  return(output)
}