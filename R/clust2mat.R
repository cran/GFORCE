# Convert K-means solution to SDP solution

#' Convert a clustering or grouping to partnership matrix.
#' 
#' Takes a clustering \eqn{G} and constructs \eqn{B(G)}.
#'
#' @param clusters length \eqn{d} vector. Assigns each variable or data point to a cluster. Cluster names
#' can be numbers or strings.
#' @return a \eqn{d x d} numeric array that contains the partnership matrix corresponding to \code{clusters}.
#'
#' @examples
#' clusters <- c(1,1,1,2,2,2,3,3)
#' B_clust <- gforce.clust2mat(clusters)
#' @export
gforce.clust2mat <- function(clusters) {
  d <- length(clusters)
  group_names <- unique(clusters)
  K <- length(group_names)
  Bhat <- matrix(rep(0,d^2),ncol=d)
  
  for(k in 1:K){
    idx <- which(clusters == group_names[k])
    group_size <- length(idx)
    for(a in idx){
      for(b in idx){
        Bhat[a,b] <- 1/group_size
      }
    }
  }
  return(Bhat)
}