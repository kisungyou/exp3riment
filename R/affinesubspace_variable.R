#' extract local information : variable dimension
#'
#' use local PCA 
#' 
#' 
#' @examples 
#' \dontrun{
#' X    = as.matrix(iris[,1:4])
#' digs = diglocal_fixed(X)
#' }
#' 
#' @export
diglocal_variable <- function(data, k=10, pcaratio=0.9, center.by.data=TRUE){
  # parameters
  mynbd  = max(1, round(k))
  myndim = max(round(ndim), 1)
  
  n = base::nrow(data)
  p = base::ncol(data)
  
  # nearest neighbor search
  nnbd = RANN::nn2(data, k=(mynbd+1))$nn.idx[,2:(mynbd+1)]
  
  # extract locality information
  list_center  = list()
  list_stiefel = list()
  
  # for (i in 1:n){
  #   # select 
  #   tgt_id  = as.vector(nnbd[i,])
  #   tgt_vec = as.vector(data[i,])
  #   tgt_mat = as.matrix(data[tgt_id,])
  #   
  #   if (center.by.data){
  #     # centering by data 
  #     tgt_centered = sweep(tgt_mat, 2, tgt_vec, FUN="-")
  #     tgt_cov      = t(tgt_centered)%*%tgt_centered
  #     
  #     list_center[[i]]  = tgt_vec
  #     list_stiefel[[i]] = rARPACK::eigs_sym(tgt_cov, myndim)$vectors
  #   } else { 
  #     # or true local PCA
  #     tgt_cov = stats::cov(tgt_mat)
  #     
  #     list_center[[i]]  = as.vector(colMeans(tgt_mat))
  #     list_stiefel[[i]] = rARPACK::eigs_sym(tgt_cov, myndim)$vectors
  #   }
  # }
  
  # return
  return(list(center=list_center, stiefel=list_stiefel))
}