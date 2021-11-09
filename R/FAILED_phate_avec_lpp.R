#' PHATE + Potential : LPP
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
#' distEB          = stats::as.dist(fast_dist(sub_data))
#' phate_native    = phate_original(distEB, n_landmark=200)$embedding
#' lpp_nbd10 = phate_avec_lpp(distEB, n_landmark=200, nbdk=10)
#' lpp_nbd20 = phate_avec_lpp(distEB, n_landmark=200, nbdk=20)
#' lpp_nbd50 = phate_avec_lpp(distEB, n_landmark=200, nbdk=50)
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,2))
#' plot(lpp_nbd10$embedding, pch=19, col=sub_label, main="LPP10-embedding")
#' plot(lpp_nbd10$projection, pch=19, col=sub_label, main="LPP10-projection")
#' 
#' plot(lpp_nbd20$embedding, pch=19, col=sub_label, main="LPP20-embedding")
#' plot(lpp_nbd20$projection, pch=19, col=sub_label, main="LPP20-projection")
#' 
#' plot(lpp_nbd50$embedding, pch=19, col=sub_label, main="LPP50-embedding")
#' plot(lpp_nbd50$projection, pch=19, col=sub_label, main="LPP50-projection")
#' par(opar)
#' }
#' 
#' @keywords internal
#' @noRd
phate_avec_lpp <- function(data, ndim=2, nbdk=5, alpha=2.0, potential=c("log","sqrt","none"), n_landmark=100){
  # --------------------- ORIGINAL PART ----------------------------------------
  # inputs
  step1_start = Sys.time()
  if (inherits(data, "dist")){
    D = data
    N = round((sqrt(8*length(D)+1)+1)/2)
  } else if (is.matrix(data)){
    N = base::nrow(data)
    D = aux_dist(data)
  } else {
    stop("* phate_avec_lpp : 'dist' or 'matrix', one of two permitted.")
  }
  M = round(n_landmark) # compression down to (M x M)
  opt_potential = match.arg(potential)   # TUNABLE PARAMETER 
  
  fun_cmds = utils::getFromNamespace("hidden_cmds", "maotai")
  fun_mmds = utils::getFromNamespace("hidden_mmds", "maotai")
  
  time_step1 = Sys.time()-step1_start
  print(paste0("* _avec_lpp : ",round(as.double(time_step1),5), "secs : step 1 complete."))
  
  # build kernel matrix & markov kernel ----------------------------------------
  step2_start = Sys.time()
  
  #mat_kernel = aux_kernel_standard(D, round(nbdk), as.double(alpha))
  #mat_P      = base::diag(1/base::rowSums(mat_kernel))%*%mat_kernel 
  mats_run   = src_standard_kernel(as.matrix(D), round(nbdk), as.double(alpha))
  mat_kernel = mats_run$mat_kernel
  mat_P      = mats_run$mat_markov
  
  time_step2 = Sys.time()-step2_start
  print(paste0("* _avec_lpp : ",round(as.double(time_step2),5), "secs : step 2 complete."))
  
  
  if (N > M){  # reduced version
    start_landmark1 = Sys.time()
    # clustering
    pseudo_data = mat_P%*%matrix(stats::rnorm(N*50), ncol=50)
    clustering  = as.integer(factor(stats::kmeans(pseudo_data, M)$cluster))
    
    # clustering  = cluster::pam(D, M)$clustering
    
    time_landmark1 = Sys.time()-start_landmark1
    print(paste0("* _avec_lpp : ",round(as.double(time_landmark1),5), "secs : landmark : (1) clustering.."))
    
    
    # compression
    start_landmark2 = Sys.time()
    compression = aux_compression(mat_P, mat_kernel, M, clustering)
    P_MN = compression$P_MN
    P_NM = compression$P_NM
    P_MM = P_MN%*%P_NM
    
    time_landmark2 = Sys.time()-start_landmark2
    print(paste0("* _avec_lpp : ",round(as.double(time_landmark2),5), "secs : landmark : (2) compression.."))
    
    
    # optimal t
    start_landmark3 = Sys.time()
    opt.t = aux_entropyrule_markov(P_MM)
    
    time_landmark3 = Sys.time()-start_landmark3
    print(paste0("* _avec_lpp : ",round(as.double(time_landmark3),5), "secs : landmark : (3) entropy decision.."))
    
    # embedding
    start_landmark4 = Sys.time()
    Pout  = P_MM
    for (i in 1:(opt.t-1)){
      Pout = P_MM%*%Pout
    }
    Pout = P_NM%*%Pout
  } else {     # naive version
    # optimal t
    opt.t = aux_entropyrule(mat_kernel)
    print(paste0("optimal t=",opt.t))
    
    # embedding
    Pout  = mat_P
    for (i in 1:(opt.t-1)){
      Pout = mat_P%*%Pout
    }
  }

  # ---------------------------- LPP Part --------------------------------------
  # Potential Coordinates
  if (all(opt_potential=="log")){
    pseudo_X = base::log(Pout+(1e-8))
  } else if (all(opt_potential=="sqrt")){
    pseudo_X = base::sqrt(Pout)
  } else {
    pseudo_X = Pout
  }
  pseudo_dist = cpp_distance(pseudo_X)
  
  # construct the affinity matrix 
  N = base::nrow(pseudo_dist)
  matW = array(0,c(N,N))
  for (n in 1:N){
    tgt_vec = as.vector(pseudo_dist[n,])
    tgt_idx = order(tgt_vec)[2:(round(nbdk)+1)]
    matW[n,tgt_idx] = 1
    matW[tgt_idx,n] = 1
  }
  matD = diag(rowSums(matW))
  
  termL   = t(pseudo_X)%*%(matD-matW)%*%pseudo_X
  termR   = t(pseudo_X)%*%matD%*%pseudo_X
  geigLR  = geigen::geigen(termL, termR)
  ntarget = length(geigLR$values)
  projmat = geigLR$vectors[,ntarget:(ntarget-round(ndim)+1)]
  projqr  = qr.Q(qr(projmat))
  
  # Return?
  output = list()
  output$embedding  = pseudo_X%*%projmat
  output$projection = pseudo_X%*%projqr
  return(output)
}