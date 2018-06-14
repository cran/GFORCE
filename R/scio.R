# estimate a columnwise SCIO estimator like \hat\Theta_{\cdot K}
# see Liu and Luo (2012)


#' SCIO Estimator.
#'
#' Estimate the precision matrix with the SCIO estimator. This algorithm is due to Liu and Luo (2012).
#' The implementation follows the active set strategy also used in the SCIO package.
#' 
#' @param C a \eqn{d x d} numeric matrix. This is the matrix of which we seek the inverse.
#' @param lambda a numeric. This is the sparsity penalty parameter.
#' @param k an integer. Indicates the column of the inverse to compute.
#' @param eps a numeric. A threshold used as a stopping criterion.
#' @param max_iter an integer. The max number of iterations of the SCIO algorithm.
#' @param R_only logical expression. If \code{R_only == FALSE}, then the included
#' native code implementation will be used. Otherwise, an R implementation is used.
#' @return a \eqn{d} dimensional numeric vector that is the \code{k}th column of the inverse of \code{C}.
#'
#' @examples
#' C <- diag(5)
#' theta_1 <- gforce.scio(C,0.01,1)
#'
#' @references T. Cai, W. Liu and X. Luo. A constrained l1 minimization approach to sparse precision
#'             matrix estimation. \emph{Journal of the American Statistical Association}, 2011.
#' 
#' @useDynLib GFORCE scio_column_R
#' @export
gforce.scio <- function(C, lambda, k = NULL, eps = 10^-6, max_iter = 10000,R_only=FALSE) {
  res <- NULL
  if(!is.null(k)) {
    if(R_only) {
      res <- scio_column(C,k,lambda,eps,max_iter)
    } else {
      d <- ncol(C)
      C_res <- .C(scio_column_R,
            C = as.double(C),
            d = as.integer(d),
            k = as.integer(k - 1),
            theta_k = numeric(d),
            lamba = as.double(lambda),
            eps = as.double(eps),
            max_iter = as.integer(max_iter))
      res <- C_res$theta_k
    }
  }

  return(res)
}

scio_column <- function(C, k, lambda, eps = 10^-6, max_iter = 10000) {
  # soft thresholding
  soft_threshold <- function(x){
    #find soft threshold
    x_th <- abs(x)-lambda

    if(x_th >= 0){
      return(sign(x)*x_th)
    }
    return(0)
  }
  d <- nrow(C)

  eps50 <- 50*eps
  done_eps <- FALSE
  theta_t <- rep(0,d)
  #warm start for theta
  theta_t[k] = C[k,k]
  nCtheta <- -C%*%theta_t

  active_set_mode <- FALSE
  active_idx <- rep(TRUE,d)

  niter <- 1
  while(niter < max_iter && !done_eps){
    delta_t <- 0
    for(i in 1:d) {
      if(!active_set_mode || active_idx[i]) {
        theta_t_i <- theta_t[i]
        theta_tp1_i <- nCtheta[i] + C[i,i]*theta_t_i
        if(i == k){
          theta_tp1_i <- theta_tp1_i + 1
        }

        #perform the update
        theta_t[i] <- soft_threshold(theta_tp1_i) / C[i,i]

        # if theta_t[i] changes, we need to update deltas and nCtheta
        if(theta_t_i != theta_t[i]){
          delta <- theta_t[i] - theta_t_i
          delta_t <- max(delta_t,abs(delta))
          nCtheta <- nCtheta - delta*C[,i]
        }

      }
    }

    #update active sets
    if(active_set_mode) {
      if(delta_t < eps) {
        active_set_mode <- FALSE
      }
    } else {
      if(delta_t < eps) {
        done_eps <- TRUE
      }
      if(delta_t < eps50) {
        active_set_mode <- TRUE
        active_idx <- theta_t != 0.0
      }
    }

    niter <- niter+1
  }

  return(theta_t)
}

scio_objective_function <- function(Chat,theta,k){
  return((1/2)*(theta%*%Chat%*%theta) - theta[k])
}

scio_likelihood <- function(C, Cinv) {
  ot <- as.numeric(unlist(determinant(Cinv)))
  if (ot[2]<=0) warning("Precision matrix estimate is not positive definite!")
  tmp <- (sum(diag(C%*%Cinv))  - ot[1] - dim(Cinv)[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}