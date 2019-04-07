# FUNCTIONS TO CONSTRUCT OPTIMALITY CERTIFICATES FOR THE
# K-MEANS SDP

#' FORCE optimality certificate.
#'
#' Given a proposed integer solution to the \eqn{K}-means SDP relaxation, this
#' function attempts to construct a solution to the dual problem with matching
#' objective value.
#' @param sol vector of length \eqn{d}. This contains the assignment of variables or
#' points to clusters.
#' @param D \eqn{d x d} matrix.
#' @param eps1 a scalar. It controls the stopping condition for the dual solution search.
#' @param eps2 a scalar. It controls the infeasibility tolerance for the dual solution to allow for numerical imprecision.
#' @param Y_T_min a scalar. THe smallest \eqn{Y_T} that the function can return. Must be greater than 0.
#'
#' @return An object with the following components:
#' \describe{
#' \item{\code{Y_T}}{ a numeric. The value of the variable \code{Y_T} in the dual solution found.}
#' \item{\code{Y_a}}{ a \eqn{d} dimensional numeric vector. The value of the variable \code{Y_a} in the dual solution found.}
#' \item{\code{feasible}}{an integer. 1 signifies that \code{sol} is optimal, 0 otherwise.}
#' }
#'
#' @examples
#' K <- 5
#' n <- 50
#' d <- 50
#' dat <- gforce.generator(K,d,n,3,graph='scalefree')
#' sig_hat <- (1/n)*t(dat$X)%*%dat$X
#' gam_hat <- gforce.Gamma(dat$X)
#' D <- diag(gam_hat) - sig_hat
#' dual_cert <- gforce.certify(dat$group_assignments,D)
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#'
#' @useDynLib GFORCE kmeans_dual_solution_primal_min_R
#' @export
gforce.certify <- function(sol,D,eps1 = 0.01,eps2 = 10^-7,Y_T_min = 0.01) {
  # check inputs
  Y_T_min <- as.double(Y_T_min)
  eps1 <- as.double(eps1)
  eps2 <- as.double(eps2)
  if(Y_T_min <= 0){
    stop('gforce.certify -- Y_T_min must be greater than 0.')
  }
  if(eps1 <= 0){
    stop('gforce.certify -- eps1 must be greater than 0.')
  }
  if(eps2 < 0){
    stop('gforce.certify -- eps1 must be greater than or equal to 0.')
  }

  K <- length(unique(sol))
  d <- dim(D)[1]
  C_result <- .C(kmeans_dual_solution_primal_min_R,
               sol=as.integer(sol),
               D=as.double(D),
               K= as.integer(K),
               dimension=as.integer(d),
               eps1 = eps1,
               eps2 = eps2,
               Y_T_min = Y_T_min,
               Y_a = numeric(d),
               Y_T = as.double(0),
               feasible = as.integer(0))
  res <- NULL
  res$Y_T <- C_result$Y_T
  res$Y_a <- C_result$Y_a
  res$feasible <- C_result$feasible
  return(res)
}


#' FORCE optimality certificate (\eqn{K} is unknown).
#'
#' Given a proposed integer solution to the adaptive \eqn{K}-means SDP relaxation, this
#' function attempts to construct a solution to the dual problem with matching
#' objective value.
#' @param sol vector of length \eqn{d}. This contains the assignment of variables or
#' points to clusters.
#' @param D \eqn{d x d} matrix.
#' @param eps1 a scalar. It controls the infeasibility tolerance for the dual solution to allow for numerical imprecision.
#'
#' @return An object with the following components:
#' \describe{
#' \item{\code{Y_a}}{ a \eqn{d} dimensional numeric vector. The value of the variable \code{Y_a} in the dual solution found.}
#' \item{\code{feasible}}{an integer. 1 signifies that \code{sol} is optimal, 0 otherwise.}
#' }
#'
#' @examples
#' K <- 5
#' n <- 50
#' d <- 50
#' dat <- gforce.generator(K,d,n,3,graph='scalefree')
#' sig_hat <- (1/n)*t(dat$X)%*%dat$X
#' gam_hat <- gforce.Gamma(dat$X)
#' D <- diag(gam_hat) - sig_hat
#' dual_cert <- gforce.certify_adapt(dat$group_assignments,D)
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#'
#' @useDynLib GFORCE kmeans_dual_solution_primal_min_nok_R
#' @export
gforce.certify_adapt <- function(sol,D,eps1 = 10^-7) {
  # check inputs
  eps1 <- as.double(eps1)

  if(eps1 <= 0){
    stop('gforce.certify_adapt -- eps1 must be greater than 0.')
  }

  K <- length(unique(sol))
  d <- dim(D)[1]
  C_result <- .C(kmeans_dual_solution_primal_min_nok_R,
               sol=as.integer(sol),
               D=as.double(D),
               K= as.integer(K),
               dimension=as.integer(d),
               eps1 = eps1,
               Y_a = numeric(d),
               feasible = as.integer(0))
  res <- NULL
  res$Y_a <- C_result$Y_a
  res$feasible <- C_result$feasible
  return(res)
}


# This function is used now for testing the C implementation
# construct dual solutions to kmeans SDP relaxation
# ga_hat -- proposed primal integer solution
# D -- objective value function in **maximization** version of the relaxation
dual_solution <- function(ga_hat,D,eps = 0.01,eps2 = 10^-7,Y_T_min = 0.01){

  return_dual <- NULL
  group_ids <- unique(ga_hat)
  K <- length(group_ids)
  d <- dim(D)[1]
  group_idxs <- list()
  group_sizes <- rep(0,K)
  group_sums <- rep(0,K)
  for(i in 1:K){
    ga_hat_k <- which(ga_hat == group_ids[i])
    group_idxs[[i]] <- ga_hat_k
    group_sizes[i] <- length(ga_hat_k)
    group_sums[i] <- sum(D[ga_hat_k,ga_hat_k])
  }

  # get primal optimal
  B_opt <- gforce.clust2mat(ga_hat)
  primal_value <- sum(D*B_opt)

  # construct dual optimal solution by binary search
  # on Y_T
  M <- matrix(rep(0,d*(d+1)),ncol=d)
  b_base <- matrix(rep(0,d+1),ncol=1)

  for(k in 1:K){
    g_k_idx <- group_idxs[[k]]
    g_k_size <- group_sizes[k]
    M[g_k_idx,g_k_idx] <- matrix(rep(1,g_k_size^2),ncol=g_k_size) + g_k_size*diag(g_k_size)
    if(g_k_size == 1){
      b_base[g_k_idx] <- D[g_k_idx,g_k_idx]
    } else{
      b_base[g_k_idx] <- colSums(D[g_k_idx,g_k_idx])
    }
  }
  M[d+1,1:d] <- 2*rep(1,d)
  b_base[d+1] <- primal_value

  Y_T_max <- abs(primal_value)
  dual_best <- NULL
  found_feasible <- FALSE

  # perform the search
  while(Y_T_max / Y_T_min > 1+eps){
    # construct b vector for this iteration
    Y_T_new <- (Y_T_max + Y_T_min) / 2
    b_new <- rep(0,d+1)
    b_new[1:d] <- b_base[1:d] - Y_T_new
    b_new[d+1] <- b_base[d+1] - Y_T_new*K

    # solve linear system
    Y_a_new <- qr.solve(M,b_new)

    # update current solution
    A <- Y_T_new*diag(d)
    for(a in 1:d){
      A <- A + Y_a_new[a]*R_a(d,a)
    }
    Y_ab_new <- A - D
    Y_ab_new[B_opt > 0] <- 0
    R_new <- A - D - Y_ab_new
    R_evals <- eigen(R_new,only.values=TRUE)$values
    em_new <- min(Re(R_evals))

    # update bounds on search
    if(em_new > -eps2){
      found_feasible <- TRUE
      dual_best$Y_T <- Y_T_new
      dual_best$Y_a <- Y_a_new
      dual_best$Y_ab <- Y_ab_new
      dual_best$R <- R_new
      Y_T_max <- Y_T_new
    } else {
      Y_T_min <- Y_T_new
    }
  }

  # clean up final solution if needed to ensure feasibility
  if(found_feasible){
    return_dual$Y_a <- dual_best$Y_a
    return_dual$Y_ab <- dual_best$Y_ab
    return_dual$Y_ab[return_dual$Y_ab < 0] = 0
    A <- dual_best$Y_T*diag(d)
    for(a in 1:d){
      A <- A + return_dual$Y_a[a]*R_a(d,a)
    }
    R <- A - D - return_dual$Y_ab
    R_evals <- eigen(R,only.values=TRUE)$values
    return_dual$Y_T <- dual_best$Y_T - min(0,min(Re(R_evals)))
    return_dual$R <- R + (return_dual$Y_T - dual_best$Y_T)*diag(d)
    return_dual$primal_value <- primal_value
    return_dual$dual_value <- 2*sum(return_dual$Y_a) + K*return_dual$Y_T
  }

  return(return_dual)
}


certify_adapt <- function(ga_hat,D,eps2 = 10^-7){
  group_ids <- unique(ga_hat)
  K <- length(group_ids)
  d <- dim(D)[1]
  group_idxs <- list()
  group_sizes <- rep(0,K)
  group_sums <- rep(0,K)
  for(i in 1:K){
    ga_hat_k <- which(ga_hat == group_ids[i])
    group_idxs[[i]] <- ga_hat_k
    group_sizes[i] <- length(ga_hat_k)
    group_sums[i] <- sum(D[ga_hat_k,ga_hat_k])
  }

  # get primal optimal
  B_opt <- gforce.clust2mat(ga_hat)
  primal_value <- sum(-D*B_opt) # because max_B <-D,B>

  M <- matrix(rep(0,d*(d+1)),ncol=d)
  b_base <- matrix(rep(0,d+1),ncol=1)

  for(k in 1:K){
    g_k_idx <- group_idxs[[k]]
    g_k_size <- group_sizes[k]
    M[g_k_idx,g_k_idx] <- matrix(rep(1,g_k_size^2),ncol=g_k_size) + g_k_size*diag(g_k_size)
    if(g_k_size == 1){
      b_base[g_k_idx] <- -D[g_k_idx,g_k_idx]
    } else{
      b_base[g_k_idx] <- colSums(-D[g_k_idx,g_k_idx])
    }
  }
  M[d+1,1:d] <- 2*rep(1,d)
  b_base[d+1] <- primal_value
  Y_a_new <- qr.solve(M,b_base)

  # update current solution
  A <- matrix(rep(0,d*d),ncol=d)
  for(a in 1:d){
    A <- A + Y_a_new[a]*R_a(d,a)
  }
  Y_ab_new <- A + D
  Y_ab_new[B_opt > 0] <- 0
  R_new <- A + D - Y_ab_new
  R_evals <- eigen(R_new,only.values=TRUE)$values
  em_new <- min(Re(R_evals))

  res <- NULL
  res$Y_ab <- Y_ab_new
  res$Y_a <- Y_a_new
  res$R_evals <- R_evals
  res$R_min <- em_new
  res$feasible <- (em_new > -eps2) && (min(Y_ab_new) > -eps2)

  return(res)
}
