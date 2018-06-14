
#' Confidence Intervals for Estimation in G-Latent Models.
#'
#' Estimate the precision matrix and construct confidence intervals for the latent and group averages graphs. The user
#' can either provide an estimate of the relevant covariance matrix or the data and cluster assignments. If cross validation
#' is selected, the data and clusters must be provided.
#' 
#' @param C_hat a \eqn{d x d} matrix. This is the estimated latent or group averages covariance matrix.
#' @param X_vals a \eqn{n x d} matrix. This is a matrix where each row is a sample from the model.
#' @param clusters a \eqn{d} dimensional integer vector. Contains the assignment of variables to groups.
#' @param alpha a numeric. \code{alpha} is the confidence level.
#' @param graph a string. It can either have value \code{'latent'} or \code{'averages'}.
#' @param variance_estimator a string. It can either have value \code{'simple'} or \code{'exact'}.
#' @param use_cv logical expression. Indicates whether or not to use cross validation.
#' @param cv_opts an object. Contains options for cross validation procedure.
#' @param lambda1 a numeric. Parameter for the first optimization problem.
#' @param lambda2 a numeric. Parameter for the second optimization problem.
#' 
#'
#' @return a \eqn{d x d x 3} array. The first \eqn{d x d} slice is the lower confidence bound, the 
#' second the point estimate, and the third the upper confidence bound.
#'
#' @examples
#' K <- 5
#' n <- 50
#' d <- 50
#' dat <- gforce.generator(K,d,n,3,graph='scalefree')
#' th_tilde <- gforce.glatent_confints(X_vals = dat$X,clusters = dat$group_assignments,
#'                                     use_cv = TRUE,graph='latent')
#'
#' @references C. Eisenach, F. Bunea, Y. Ning, and C. Dinicu. Efficient, High-Dimensional Inference for Cluster-Based Graphical Models. \emph{Manuscript submitted for publication}, 2018. 
#'
#' @seealso \code{\link{gforce.glatent_confints.cv_defaults}}
#' @export
gforce.glatent_confints <- function(C_hat=NULL,X_vals=NULL,clusters=NULL,alpha = 0.05,
                                    graph='latent',variance_estimator='simple',use_cv = FALSE,
                                    cv_opts = NULL,lambda1=NULL,lambda2 = NULL) {
  res <- NULL
  # if user specified C_hat
  if(!is.null(C_hat)){
    if(graph == 'latent') {
      # cross validation cannot be used
      stop('gforce.glatent_confints -- Option not implemented.')
    } else if(graph == 'averages') {
      # cross validation cannot be used
      stop('gforce.glatent_confints -- Option not implemented.')
    } else {
      stop('gforce.glatent_confints -- You must specify either latent or averages graphs.')
    }

  } else if(!is.null(clusters) && !is.null(X_vals)) {
    if(graph == 'latent') {
      if(use_cv){
          res <- latent_confidence_intervals_all_cv(X_vals,clusters,alpha,cv_opts,variance_estimator=variance_estimator)
      } else {
        stop('gforce.glatent_confints -- Option not implemented.')
      }
    } else if(graph == 'averages') {
      if(use_cv){
        res <- averages_confidence_intervals_all_cv(X_vals,clusters,alpha,cv_opts)
      } else {
        stop('gforce.glatent_confints -- Option not implemented.')
      }
    } else {
      stop('gforce.glatent_confints -- You must specify either latent or averages graphs.')
    }

  } else {
    stop('gforce.glatent_confints -- You must specify one of C_hat or clusters and X.')
  }
  return(res)
}

#' Default Cross Validation Options for Confidence Intervals.
#'
#' @return An object with following components
#' \describe{
#' \item{\code{grid_density}}{an integer. Indicates the number of values in the search grid.}
#' \item{\code{lambda1_min_exp}}{a numeric. Minimum exponent for values of \eqn{\lambda_1} to search over.}
#' \item{\code{lambda1_max_exp}}{a numeric. Maximum exponent for values of \eqn{\lambda_1} to search over.}
#' \item{\code{lambda2_min_exp}}{a numeric. Minimum exponent for values of \eqn{\lambda_2} to search over.}
#' \item{\code{lambda2_max_exp}}{a numeric. Maximum exponent for values of \eqn{\lambda_2} to search over.}
#' \item{\code{num_folds}}{an integer. The number of cross-validation folds to use.}
#' }
#' @seealso \code{\link{gforce.glatent_confints}}
#' @export
gforce.glatent_confints.cv_defaults <- function() {
  res <- NULL
  res$grid_density <- 5
  res$lambda1_min_exp <- -8
  res$lambda1_max_exp <- 0
  res$lambda2_min_exp <- -8
  res$lambda2_max_exp <- 0
  res$num_folds <- 5
  return(res)
}

latent_confidence_intervals_all_cv <- function(X_vals,group_assignments,alpha,cv_opts=NULL,variance_estimator='simple') {
  # check for cv_opts
  if(is.null(cv_opts)){
    cv_opts <- gforce.glatent_confints.cv_defaults()
  }

  # ensure group labels ordered properly
  group_assignments <- order_group_assignments(group_assignments)
  
  # create C hat estimator
  n <- nrow(X_vals)
  groups <- unique(group_assignments)
  K <- length(groups)
  res <- array(0,c(K,K,3))
  sig_hat <- (t(X_vals)%*%X_vals)/n
  gam_hat <- glatent_Gamma_hat(sig_hat,group_assignments)
  Chat <- latent_spectrum_check(C_hat(sig_hat,gam_hat,group_assignments))
  
  #estimate theta *columnwise*
  theta_hat <- NULL
  theta_hat <- array(0,c(K,K))
  for(k in 1:K){
    objective_function <- function(ch,tha) scio_objective_function(ch,tha,k)
    function_solver <- function(ch,lambda) gforce.scio(ch,lambda,k)
    lambda1 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,
                                    cv_opts$lambda1_min_exp,cv_opts$lambda1_max_exp,cv_opts$num_folds,graph='latent')
    theta_hat[,k] <- gforce.scio(Chat,lambda1,k)
  }

  #estimate v *row-wise*
  v_hat <- array(0,c(K,K))
  for(t in 1:K){
    objective_function <- function(ch,tha) max(abs(ch%*%tha))
    function_solver <- function(ch,lambda) v_hat(ch,t,lambda)
    lambda2 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,
                                    cv_opts$lambda2_min_exp,cv_opts$lambda2_max_exp,cv_opts$num_folds,graph='latent')
    v_hat[t,] <- v_hat(Chat,t,lambda2)
  }
  #estimate \tilde \theta_t,k *row-wise*
  for(t in 1:K){
    res[t,,2] <- decorrelated_estimator_t(Chat,theta_hat,v_hat[t,],t)
  }
  
  #construct upper and lower confidence bounds
  z_score <- stats::qnorm(1 - (alpha/2))

  if(variance_estimator=='simple') {
    for(t in 1:K){
      for(k in 1:K){
        variance <- (theta_hat[t,k]^2 + theta_hat[t,t]*theta_hat[k,k])#/(xi_hat[t,t]^2)
        conf_range <- sqrt(variance)*z_score/sqrt(n)
        res[t,k,1] <- res[t,k,2] - conf_range
        res[t,k,3] <- res[t,k,2] + conf_range
      }
    }
  } else if(variance_estimator=='exact') {
    latent_variances <- latent_graph_variances(theta_hat,v_hat,gam_hat,group_assignments)
    for(k in 1:K){
      for(t in 1:K){
        conf_range <- sqrt(latent_variances[t,k])*z_score/sqrt(n)
        res[t,k,1] <- res[t,k,2] - conf_range
        res[t,k,3] <- res[t,k,2] + conf_range
      }
    }
  } else {
    stop('Variance Estimator must be simple or exact')
  }

  return(res)
}

averages_confidence_intervals_all_cv <- function(X_vals,group_assignments,alpha,cv_opts=NULL) {
  # check for cv_opts
  if(is.null(cv_opts)){
    cv_opts <- gforce.glatent_confints.cv_defaults()
  }

  # ensure group labels ordered properly
  group_assignments <- order_group_assignments(group_assignments)
  
  # create C hat estimator
  n <- nrow(X_vals)
  groups <- unique(group_assignments)
  K <- length(groups)
  res <- array(0,c(K,K,3))
  sig_hat <- (t(X_vals)%*%X_vals)/n
  s_hat <- S_hat(sig_hat,group_assignments)
  
  #estimate theta *columnwise*
  xi_hat <- array(0,c(K,K))
  for(k in 1:K){
    objective_function <- function(ch,tha) scio_objective_function(ch,tha,k)
    function_solver <- function(ch,lambda) gforce.scio(ch,lambda,k)
    lambda1 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,
                                    cv_opts$lambda1_min_exp,cv_opts$lambda1_max_exp,cv_opts$num_folds,graph='averages')
    xi_hat[,k] <- gforce.scio(s_hat,lambda1,k)
  }
  
  #estimate v *row-wise*
  v_hat <- array(0,c(K,K))
  for(t in 1:K){
    objective_function <- function(ch,tha) max(abs(ch%*%tha))
    function_solver <- function(ch,lambda) v_hat(ch,t,lambda)
    lambda2 <- cv_lambda_selection(objective_function,function_solver,X_vals,group_assignments,
                                    cv_opts$lambda2_min_exp,cv_opts$lambda2_max_exp,cv_opts$num_folds,graph='averages')
    v_hat[t,] <- v_hat(s_hat,t,lambda2)
  }
  #estimate \tilde \theta_t,k *row-wise*
  for(t in 1:K){
    res[t,,2] <- decorrelated_estimator_t(s_hat,xi_hat,v_hat[t,],t)
  }
  
  #construct upper and lower confidence bounds
  z_score <- stats::qnorm(1 - (alpha/2))
  for(t in 1:K){
    for(k in 1:K){
      variance <- (xi_hat[t,k]^2 + xi_hat[t,t]*xi_hat[k,k])#/(xi_hat[t,t]^2)
      conf_range <- sqrt(variance)*z_score/sqrt(n)
      res[t,k,1] <- res[t,k,2] - conf_range
      res[t,k,3] <- res[t,k,2] + conf_range
    }
  }
  
  return(res)
}

