#' Graphical PHATE for Effective Resistance
#' 
#' 
#' 
#' 
#' 
#' 
#' @export
ver2_phate_graphE <- function(input, ndim=2, nbdk=5, alpha=2.0, time_step=NULL, add_diagonal=TRUE){
  # parameters
  opt_algorithm = "mmds"
  opt_potential = "log"
  if ((is.null(time_step))&&(length(time_step)<1)){
    markov_rule = TRUE # use entropy-based determination
    markov_step = 1
  } else {
    markov_rule = FALSE
    markov_step = max(1, round(time_step))
  }
  
  # STEP 1. Use Adjacency Matrix Directly
  mat_input = aux_binarynetwork(input)
  if (add_diagonal){
    diag(mat_input) = 1
  }
  if (isSymmetric(mat_input)){
    ER = aux_effectivesym(mat_input)
  } else {
    ER = aux_effective(mat_input)
  }
  ER[(ER<=0)] = 0
  D = stats::as.dist(base::sqrt(ER))
  
  # STEP 2. Pass onto 'ver2_phate_original' to do the rest
  output = ver2_phate_original(D, ndim=ndim, nbdk=nbdk, alpha=alpha, time_step=time_step)
  return(output)
}

# Run some examples
# # load the data
# data("dolphins", package="exp3riment")
# dat_graph = dolphins$igraph
# dat_label = dolphins$label
# 
# dol_nodiag = ver2_phate_graphE(dat_graph, add_diagonal = FALSE)
# dol_oodiag = ver2_phate_graphE(dat_graph, add_diagonal = TRUE)
# 
# require(igraph)
# par(mfrow=c(1,3))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA, 
#      main="igraph")
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA, 
#      main="PHATE+nodiag", layout=dol_nodiag$embedding)
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA, 
#      main="PHATE+oodiag", layout=dol_oodiag$embedding)