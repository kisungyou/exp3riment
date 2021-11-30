#' Distance using Conway Embedding
#' 
#' 
#' 
#' @examples 
#' \dontrun{
#' # example with Iris data
#' X    = as.matrix(iris[,1:4])
#' lab  = as.factor(iris[,5])
#' digs = diglocal_fixed(X, center.by.data=FALSE)
#' 
#' # compute the distance
#' dobj = dist_conway(digs$center, digs$stiefel)
#' 
#' # embed by MDS
#' embed2_data   = cmdscale(stats::dist(X), k=2)
#' embed2_conway = cmdscale(dobj, k=2)
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
dist_conway <- function(list_center, list_stiefel){
  # input controls
  if (!is.list(list_center)){
    stop("* dist_conway : 'list_center' is not a list.")
  }
  if (!is.list(list_stiefel)){
    stop("* dist_conway : 'list_stiefel' is not a list.")
  }
  if (length(list_center)!=length(list_stiefel)){
    stop("* dist_conway : no matching distances.")
  }
  
  # compute distances
  n = length(list_center)
  p = length(list_center[[1]])
  
  projs = array(0,c(p,p,n))
  for (i in 1:n){
    tmpst = list_stiefel[[i]]
    if (is.vector(tmpst)){
      tmpst = as.matrix(tmpst)
    }
    projs[,,i] = tmpst%*%t(tmpst)
  }
  
  Dmat = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d1 <- sqrt(sum((as.vector(list_center[[i]]) - as.vector(list_center[[j]]))^2))
      d2 <- base::norm(as.matrix(projs[,,i])-as.matrix(projs[,,j]), type="F")/sqrt(2)
      Dmat[i,j] <- Dmat[j,i] <- sqrt((d1^2)+(d2^2))
    }
  }
  return(stats::as.dist(Dmat))
}