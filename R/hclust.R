# Hierarchical Clustering with cluster estimation.

#' Hierarchical Clustering with Estimation of \eqn{K}.
#'
#' Clusters \eqn{n} points of dimension \eqn{m} using a complete linkage algorithm and estimates \eqn{K}.
#' 
#' @param X \eqn{n x m} matrix. Each row is treated as a point in \eqn{R^m}.
#' @param dists \eqn{n x n} symmetric matrix. This encodes the distances between the \eqn{n} points.
#' @param R_only logical expression. If \code{R_only == FALSE}, then the included
#' native code implementation will be used. Otherwise, an R implementation is used.
#' @return Returns an object with the components:
#' \describe{
#' \item{\code{K}}{an estimate of the number of clusters.}
#' \item{\code{clusters}}{a \eqn{n} dimensional integer vector. Entry \eqn{i} to the cluster assignment of the data point given by row \eqn{i} of \code{X}.}
#' \item{\code{MSE}}{a \eqn{n} dimensional vector of the mean squared errors of each choice of \eqn{K}.}
#' }
#'
#' @examples
#' m <- 10 
#' n <- 10
#' X <- matrix(MASS::mvrnorm(m*n,rep(0,m*n),diag(m*n)), nrow = n)
#' hc_res <- gforce.hclust(X=X)
#'
#' @references D. Defays. An efficient algorithm for a complete link method. \emph{The Computer Journal}, 1977.
#'
#' @useDynLib GFORCE hclust_R
#' @export
gforce.hclust <- function(X=NULL,dists=NULL,R_only=FALSE) {
  res <- NULL
  if(R_only) {
    if(is.null(X)) {
      stop('gforce.hclust -- need to specify X in R_only mode')
    }
    dX <- stats::dist(X)
    hc <- stats::hclust(dX)
    d <- ncol(X)

    MSEs <- rep(0,d)
    for(k in 2:d){
      cc <- stats::cutree(hc,k=k)
      cc_mat <- gforce.clust2mat(cc)
      MSEs[k] <- 0.5*sum(cc_mat*as.matrix(dX))
    }

    lc <- L_curve_criterion(MSEs[2:(d-1)])
    res$K <- lc$max + 1
    res$clusters <- stats::cutree(hc,k=res$K)
    res$MSE <- MSEs
  } else {
    if(is.null(X) && is.null(dists)) {
      stop('gforce.hclust -- need to specify at least one of X and dists')
    }
    if(is.null(dists)) {
      n <- nrow(X)
      dists <- matrix(0,ncol=n,nrow=n)
      ips <- X %*% t(X)
      o <- rep(1,n)
      ips_d <- diag(ips)
      dists <- ips_d%*%t(o) + o%*%t(ips_d) - 2*ips
      dists <- sqrt(dists)
    }

    n <- nrow(dists)

    C_result <- .C(hclust_R,
            D = as.double(dists),
            n = as.integer(n),
            clusters = as.integer(rep(1,n)),
            K = as.integer(1),
            MSE = numeric(n))

    res$K <- C_result$K
    res$clusters <- C_result$clusters
    res$MSE <- C_result$MSE
  }
  
  return(res)
}


#' Hierarchical Clustering Agglomeration.
#'
#' @param X \eqn{n x m} matrix. Each row is treated as a point in \eqn{R^m}.
#' @param dists \eqn{n x n} symmetric matrix. This encodes the distances between the \eqn{n} points.
#'
#' @useDynLib GFORCE hclust_agglomerate_R
#' @export
gforce.hclust.agglomerate <- function(X=NULL,dists = NULL) {
  res <- NULL

  if(is.null(X) && is.null(dists)) {
    stop('gforce.hclust -- need to specify at least one of X and dists')
  }
  if(is.null(dists)) {
    n <- nrow(X)
    dists <- matrix(0,ncol=n,nrow=n)
    ips <- X %*% t(X)
    o <- rep(1,n)
    ips_d <- diag(ips)
    dists <- ips_d%*%t(o) + o%*%t(ips_d) - 2*ips
    dists <- sqrt(dists)
  }

  n <- nrow(dists)

  C_result <- .C(hclust_agglomerate_R,
          D = as.double(dists),
          n = as.integer(n),
          ag1 = as.integer(rep(1,n)),
          ag2 = as.integer(rep(1,n)),
          agdist = numeric(n))

  res$ag1 <- C_result$ag1 + 1
  res$ag2 <- C_result$ag2 + 1
  res$agdist <- C_result$agdist

  return(res)
}

#' Hierarchical Clustering -- Convert Agglomeration to clustering.
#' 
#' @param hc an object. This encodes the order of agglomeration.
#' @param K an integer. \code{K} is the number of clusters
#'
#' @export
gforce.hclust.agg2clust <- function(hc,K) {
  res <- NULL
  ag_min <- hc$ag1
  ag_max <- hc$ag2
  n <- length(ag_min)

  # K groups means that we need to perform n-K merges
  clusts <- 1:n
  for(i in 1:(n-K)) {
    idx1 <- ag_min[i]
    idx2 <- ag_max[i]

    clust_1 <- clusts[idx1]
    clust_2 <- clusts[idx2]
    clust_2idx <- which(clusts == clust_2)
    clusts[clust_2idx] <- clust_1
  }


  res$clusters <- clusts
  return(res)
}