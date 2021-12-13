#' My Version of MAGIC
#' 
#' (n x p) data : cellxgene count matrix : let's not take sparse ones yet.
#' 
#' 
#' @export
lab_magic <- function(cellxgene, knn=5, t=3, apply_sqrt=FALSE){
  # names of cell's and gene's
  name_cell = rownames(cellxgene)
  name_gene = colnames(cellxgene)

  
  # unscaled 3d output : n x p x (t+1)
  output = cpp_magic(cellxgene, round(t), round(knn), as.logical(apply_sqrt))
  
  # quantiles
  num_genes = base::ncol(output)
  col_quantiles = rep(0, num_genes)
  for (j in 1:num_genes){
    col_quantiles[j] = as.double(stats::quantile(as.vector(cellxgene[,j]), 0.99))
  }
  
  
  for (it in 1:(round(t+1))){ # for each slice
    for (j in 1:num_genes){
      tgt_col = as.vector(output[,j,it])
      output[,j,it] = tgt_col*col_quantiles[j]/as.double(max(tgt_col))
    }
  }
  return(output)
}