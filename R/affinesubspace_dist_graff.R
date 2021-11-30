#' Distance using Graff (affine grassmann by Lim)
#' 
#' @examples 
#' \dontrun{
#' # example with Iris data
#' X    = as.matrix(iris[,1:4])
#' lab  = as.factor(iris[,5])
#' digs = diglocal_fixed(X)
#' 
#' # compute the distance
#' dobj = dist_graff(digs$center, digs$stiefel)
#' 
#' # embed by MDS
#' embed2_data  = cmdscale(stats::dist(X), k=2)
#' embed2_grass = cmdscale(dobj, k=2)
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed2_data,   col=lab, pch=19, main="Data MDS")
#' plot(embed2_conway, col=lab, pch=19, main="Graff MDS")
#' par(opar)
#' }
#' 
#' @export
dist_graff <- function(list_center, list_stiefel){
  # input controls
  if (!is.list(list_center)){
    stop("* dist_graff : 'list_center' is not a list.")
  }
  if (!is.list(list_stiefel)){
    stop("* dist_graff : 'list_stiefel' is not a list.")
  }
  if (length(list_center)!=length(list_stiefel)){
    stop("* dist_graff : no matching distances.")
  }
  
  # compute distances
  n = length(list_center)
  p = length(list_center[[1]])
  
  # prepare a stiefel coordinate
  vec_stiefel = list()
  for (i in 1:n){
    Ab = qr.Q(qr(cbind(list_stiefel[[i]], list_center[[i]])))
    nb = base::ncol(Ab)
    
    A  = Ab[,1:(nb-1)]
    b0 = Ab[,nb]
    bc = sqrt(1+(sum(b0^2)))
    
    vec_stiefel[[i]] = rbind(cbind(A, b0/bc),  c(rep(0,nb-1),1/bc))
  }
  
  # compute the pairwise distances
  Dmat = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      svdrun  = base::svd(t(vec_stiefel[[i]])%*%vec_stiefel[[j]])
      pangles = base::acos(svdrun$d)
      pangles[is.na(pangles)] = 0
      Dmat[i,j] <- Dmat[j,i] <- sqrt(sum(pangles^2))
    }
  }
  return(stats::as.dist(Dmat))
}
