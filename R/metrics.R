# FUNCTIONS FOR CLUSTERING AND METRICS

#' Evaluates the correctness of a clustering solution.
#'
#' This can be used for either data or variable clustering. Cluster names
#' can be either strings or numbers. Arbitrary data types cannot be used.
#' @param true_clust length \eqn{d} vector of cluster assignments. This represents the
#' true, or reference, clustering.
#' @param est_clust length \eqn{d} vector of cluster assignments. This represents the
#' estimated clustering.
#' @param method the method used to evaluate the quality of the clustering solution \code{est_clust}.
#' The three options are \code{'purity'}, \code{'perfect'}, \code{'misclassified-points'}.
#' @return Returns a numeric that represents the value of the chosen metric on the two clusterings \code{true_clust} and \code{est_clust}.
#' @examples
#' clust1 <- c(1,1,1,2,2,2,3,3,3)
#' clust2 <- c(1,1,2,1,2,2,3,3,3)
#' gforce.metrics(clust1,clust2,method='purity')
#' @export
gforce.metrics <- function(true_clust,est_clust,method='purity') {
  if(method=='purity'){
    res <- purity_measure(true_clust,est_clust)
    return(res)
  } else if(method=='perfect'){
    res <- check_perfect_recovery(true_clust,est_clust)
    return(res)
  } else if(method=='misclassified-points'){
    res <- misclassified_points(true_clust,est_clust)
  } else {
    stop("gforce.metrics -- Need to specify either 'purity', 'perfect' or 'misclassified-points' as method.")
  }
  return(res)
}

check_perfect_recovery <- function(original,recovered){
  group_ids_orig <- unique(original)
  group_ids_recov <- unique(recovered)
  K1 <- length(group_ids_orig)
  K2 <- length(group_ids_recov)
  if(K1 != K2){
    return(FALSE)
  }
  K <- K1
  original_groups <- list()
  recovered_groups <- list()
  for(i in 1:K){
    original_groups[[i]]  <- which(original == group_ids_orig[i])
    recovered_groups[[i]]  <- which(recovered == group_ids_recov[i])
  }
  same <- TRUE
  for(i in 1:K){
    orig_group <- original_groups[[i]]
    found_match <- FALSE
    for(j in 1:K){
      # check for match in recovered group j
      recov_group <- recovered_groups[[j]]
      if(length(orig_group) == length(recov_group)){
        if(length(orig_group) == sum(orig_group %in% recov_group)){
          found_match <- TRUE
        }
      }
    }
    if(!found_match){
      same <- FALSE
    }
  }
  
  return(same)
}

purity_measure <- function(ga,ga_hat){
  group_ids_ga <- unique(ga)
  group_ids_ga_hat <- unique(ga_hat)
  K <- length(group_ids_ga_hat)
  K2 <- length(unique(group_ids_ga))
  n <- length(ga_hat)
  p <- 0
  for(k in 1:K){
    m <- 0
    for(kp in 1:K2){
      group_ga_hat_k <- which(ga_hat == group_ids_ga_hat[k])
      group_ga_kp <- which(ga == group_ids_ga[kp])
      int_size <- sum(group_ga_hat_k %in% group_ga_kp)
      m <- max(m,int_size)
    }
    p <- p + m
  }
  p <- p/n
  return(p)
}

misclassified_points <- function(ga,ga_hat){
  group_ids_ga <- unique(ga)
  group_ids_ga_hat <- unique(ga_hat)
  K <- length(group_ids_ga_hat)
  K2 <- length(unique(group_ids_ga))
  n <- length(ga_hat)
  p <- 0
  misclassified_points <- list()
  for(k in 1:K){
    m <- 0
    misclassified_points_from_grp <- c()
    group_ga_hat_k <- which(ga_hat == group_ids_ga_hat[k])
    for(kp in 1:K2){
      group_ga_kp <- which(ga == group_ids_ga[kp])
      corr_class <- group_ga_hat_k %in% group_ga_kp
      int_size <- sum(corr_class)
      if(int_size > m){
        misclassified_points_from_grp <- setdiff(group_ga_hat_k,group_ga_kp[corr_class])
        m <- int_size
      }
    }
    misclassified_points[[k]] <- misclassified_points_from_grp
  }
  return(misclassified_points)
}

kmeans_repeater <- function(sig,num_repeat,ga){
  av_purity <- 0
  percent_perfect <- 0
  K <- length(unique(ga))
  for(i in 1:num_repeat){
    cur_purity <- purity_measure(kmeanspp(sig,K),ga)
    av_purity <- av_purity + cur_purity
    if(cur_purity == 1){
      percent_perfect <- percent_perfect+1
    }
  }
  res <- NULL
  res$percent_perfect <- percent_perfect/num_repeat
  res$average <- av_purity/num_repeat
  return(res)
}

