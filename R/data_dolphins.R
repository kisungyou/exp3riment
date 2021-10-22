#' Load dolphins data
#' 
#' This data consists of social network of 62 dolphins. This small undirected network 
#' is known to consist of two communities of sizes 41 and 21 nodes each.
#' 
#' @usage data(dolphins)
#' 
#' @examples
#' \donttest{
#' # load the data
#' data(dolphins, package="T4network")
#' }
#' 
#' @format a named list containing\describe{
#' \item{igraph}{an \code{'igraph'} object.}
#' \item{name}{a length-\eqn{62} vector of node names for dolphins.}
#' \item{label}{a length-\eqn{62} vector of known community membership.}
#' }
#' 
#' @concept data
"dolphins"