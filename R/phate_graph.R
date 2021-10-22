#' Graphical PHATE for Asymmetric Graph
#' 
#' 
#' 
#' @examples 
#' \dontrun{
#' # load the data
#' data("dolphins", package="exp3riment")
#' dat_graph = dolphins$igraph
#' dat_label = dolphins$label
#' 
#' # try PHATE
#' dol_log_cmds = phate_graph(dat_graph, potential="log", alg="cmds", n_landmark=50)
#' dol_log_mmds = phate_graph(dat_graph, potential="log", alg="mmds", n_landmark=50)
#' 
#' # visualize
#' require("igraph")
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#'      main="igraph")
#' plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#'      main="PHATE+CMDS", layout=dol_log_cmds$embedding)
#' plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#'      main="PHATE+MMDS", layout=dol_log_mmds$embedding)
#' par(opar)
#' }
#' 
#' @export
phate_graph <- function(input, ndim=2, nbdk=5, alpha=2.0, alg=c("cmds","mmds"), potential=c("log","sqrt","none"), n_landmark=100){
  # Network Metrization
  input_mat = aux_binarynetwork(input)
  if (isSymmetric(input_mat)){
    ER = aux_effectivesym(input_mat)
  } else {
    ER = aux_effective(input_mat)
  }
  ER[(ER<=0)] = 0
  D = stats::as.dist(base::sqrt(ER))
  
  # Run Standard PHATE then
  myalg  = match.arg(alg)
  mypot  = match.arg(potential)
  output = phate_original(D, ndim=ndim, nbdk=nbdk,
                          alpha=alpha, alg=myalg, 
                          potential=mypot, n_landmark = n_landmark)
  return(output)
}



# # Comparison of Multiple Layouts ------------------------------------------
# # load the data
# data("dolphins", package="exp3riment")
# dat_graph = dolphins$igraph
# dat_label = dolphins$label
# 
# # try PHATE
# dol_log_mmds = phate_graph(dat_graph, potential="log", alg="mmds")
# 
# # visualize
# require("igraph")
# opar <- par(no.readonly=TRUE)
# par(mfrow=c(2,3))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="PHATE+MMDS", layout=dol_log_mmds$embedding)
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="igraph-MDS", layout=layout_with_mds(dat_graph))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="igraph-GEM", layout=layout_with_gem(dat_graph))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="igraph-KK", layout=layout_with_kk(dat_graph))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="igraph-FR", layout=layout_with_fr(dat_graph))
# plot(dat_graph, vertex.color=dat_label, vertex.label=NA,
#      main="igraph-Sugiyama", layout=layout_with_sugiyama(dat_graph)$layout)
# par(opar)
# 
# dd = aux_effective(as_adjacency_matrix(dat_graph))
# D  = stats::as.dist(sqrt(dd))
# plot(cmdscale(D))
