# Interface to C solver for K-means SDP

# FORCE solves the minimization form of the K-means SDP
# e.g. D can be a distances matrix or negated covariance
# INPUTS: Needs to have D and K. Optional force_opts
# optional initial feasible, optional E
# OUTPUTS:


#' FORCE \eqn{K}-means solver.
#'
#' Solves the Peng-Wei K-means SDP Relaxation using the FORCE algorithm.
#'
#' @param D a matrix \eqn{D} as defined above.
#' @param K number of clusters.
#' @param force_opts tuning parameters. \code{NULL} signifies defaults will be used.
#' @param D_Kmeans matrix to be used for initial integer solution. \code{NULL} signifies that \code{D} will be used.
#' @param X0 initial iterate. \code{NULL} signifies that it will be generated randomly from \code{D_Kmeans}. If supplied, \code{E} must be supplied as well.
#' @param E strictly feasible solutions. \code{NULL} signifies that it will be generated randomly. If supplied, \code{X0} must be supplied as well.
#' @param R_only logical expression. If \code{R_only == FALSE}, then the included
#' native code implementation will be used. Otherwise, an R implementation is used.
#'
#' @return An object with following components
#' \describe{
#' \item{\code{Z_T}}{Final iterate of the projected gradient descent algorithm run on the smoothed eigenvalue problem.}
#' \item{\code{B_Z_T}}{Projection of \code{Z_T} to the border of the positive semi-definite cone.}
#' \item{\code{B_Z_T_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_Z_T}.}
#' \item{\code{Z_best}}{Iterate with best objective value found during projected gradient descent on the smoothed eigenvalue problem.}
#' \item{\code{B_Z_best}}{Projection of \code{Z_best} to the border of the positive semi-definite cone.}
#' \item{\code{B_Z_best_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_Z_T}.}
#' \item{\code{km_best}}{Best clustering in terms of objective value of the SDP relaxation. This is found by running Lloyd's algorithm on the rows of \code{D_kmeans_matrix}.}
#' \item{\code{B_km}}{Partnership matrix corresponding to \code{km_best}.}
#' \item{\code{km_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_km}.}
#' \item{\code{km_best_time}}{Time elapsed (in seconds) until \code{km_best} was found.}
#' \item{\code{km_iter_best}}{Number of times a \eqn{K}-means algorithm was run before \code{km_best} was found.}
#' \item{\code{km_iter_total}}{Total number of calls to a \eqn{K}-means solver (such as Lloyd's algorithm).}
#' \item{\code{dual_certified}}{1 if a dual certificate was found, and 0 otherwise.}
#' \item{\code{dual_certified_grad_iter}}{Number of gradient updates performed before a dual certificate was found.}
#' \item{\code{dual_certified_time}}{Time elapsed (in seconds) until dual certificate was found for \code{B_km}.}
#' \item{\code{grad_iter_best}}{Gradient iteration where \code{Z_best} was computed.}
#' \item{\code{grad_iter_best_time}}{Time elapsed (in seconds) when the update \code{grad_iter_best} was performed.}
#' \item{\code{total_time}}{Total time elapsed (in seconds) during call to \code{gforce.FORCE}.}
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
#' res <- gforce.FORCE(D,K)
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#' @references J. Peng and Y. Wei. Approximating K-means-type Clustering via Semidefinite Programming. \emph{SIAM Journal on Optimization}, 2007.
#' @references J. Renegar.  Efficient first-order methods for linear programming and semidefinite programming. \emph{arXiv:1409.5832}, 2014.
#'
#' @seealso \code{\link{gforce.defaults}}
#' @useDynLib GFORCE FORCE_R
#' @export
gforce.FORCE <- function(D,K,force_opts = NULL,D_Kmeans = NULL, X0 = NULL,
                         E = NULL, R_only = FALSE) {
  d <- ncol(D)

  if(is.null(force_opts)){
      force_opts <- gforce.defaults(d)
  }

  if(is.null(D_Kmeans)) {
      D_Kmeans <- D
  }

  if(is.null(X0) && is.null(E)){
    km_res <- gforce.kmeans(-D_Kmeans,K,R_only)
    km_res <- km_res$clusters
    km_sol <- gforce.clust2mat(km_res)
    o <- rep(1,d)
    a <- (K-1)/(d-1)
    b <- (d-K)/(d^2-d)
    E <- a*diag(d) + b*o%*%t(o)
    X0 <- force_opts$initial_mixing*km_sol + (1-force_opts$initial_mixing)*E
  } else if(is.null(X0)){
    stop('FORCE -- Either specify both X0 and E or neither\r\n')
  } else if(is.null(E)){
    stop('FORCE -- Either specify both X0 and E or neither\r\n')
  }

  # make sure that initialization is well defined
  if(sum(X0*D) >= sum(E*D)) {
    stop('gforce.FORCE -- D^TX_0 >= D^TE. Check D, and if D is well specified, pass valid X0 as argument.')
  }

  E_EVEV <- eigen(E)
  E_V <- E_EVEV$vectors
  E_D <- diag(E_EVEV$values)
  E_sqrt <- E_V%*%(E_D^(0.5))%*%t(E_V)
  ESI <- solve(E_sqrt)

  res <- NULL

  if(!R_only){
    C_result <- .C(FORCE_R,
            D = as.double(D),
            D_Kmeans = as.double(D_Kmeans),
            E = as.double(E),
            ESI = as.double(ESI),
            X0 = as.double(X0),
            d = as.integer(d),
            K = as.integer(K),
            verbosity = as.integer(force_opts$verbose),
            kmeans_iter = as.integer(force_opts$kmeans_iter),
            dual_frequency = as.integer(force_opts$dual_frequency),
            max_iter = as.integer(force_opts$max_iter),
            finish_pgd = as.integer(force_opts$finish_pgd),
            primal_only = as.integer(force_opts$primal_only),
            number_restarts = as.integer(length(force_opts$restarts)),
            restarts = as.integer(force_opts$restarts),
            alpha = as.double(force_opts$alpha),
            eps_obj = as.double(force_opts$eps_obj),
            early_stop_mode = as.integer(force_opts$early_stop_mode),
            early_stop_lag = as.integer(force_opts$early_stop_lag),
            early_stop_eps = as.double(force_opts$early_stop_eps),
            Z_T = numeric(d^2),
            B_Z_T = numeric(d^2),
            Z_T_lmin = as.double(1.0),
            Z_best = numeric(d^2),
            B_Z_best = numeric(d^2),
            Z_best_lmin = as.double(2.0),
            B_Z_T_opt_val = as.double(3.0),
            B_Z_best_opt_val = as.double(4.0),
            km_opt_val = as.double(5.0),
            km_best = as.integer(numeric(d)),
            km_best_time = as.double(6.0),
            km_iter_best = as.integer(7),
            km_iter_total = as.integer(8),
            dc = as.integer(9),
            dc_time = as.double(10.0),
            dc_grad_iter = as.integer(11),
            grad_iter_best = as.integer(12),
            grad_iter_best_time = as.double(13.0),
            total_time = as.double(14))

    # Build Result
    res$Z_T <- matrix(C_result$Z_T,ncol=d)
    res$B_Z_T <- matrix(C_result$B_Z_T,ncol=d)
    res$B_Z_T_opt_val <- C_result$B_Z_T_opt_val
    res$Z_best <- matrix(C_result$Z_best,ncol=d)
    res$B_Z_best <- matrix(C_result$B_Z_best,ncol=d)
    res$B_Z_best_opt_val <- C_result$B_Z_best_opt_val
    res$km_best <- C_result$km_best
    res$B_km <- gforce.clust2mat(res$km_best)
    res$km_opt_val <- C_result$km_opt_val
    res$km_best_time <- C_result$km_best_time
    res$km_iter_best <- C_result$km_iter_best
    res$km_iter_total <- C_result$km_iter_total
    res$dual_certified <- C_result$dc
    res$dual_certified_grad_iter <- C_result$dc_grad_iter
    res$dual_certified_time <- C_result$dc_time
    res$grad_iter_best <- C_result$grad_iter_best
    res$grad_iter_best_time <- C_result$grad_iter_best_time
    res$total_time <- C_result$total_time
  } else {
    res <- primal_dual_smoothed_adar(D,D_Kmeans,K,force_opts,X0,E,ESI)
    cat("RETURNED TO FORCE")
  }
  res$X0 <- X0
  res$E <- E
  res <- res[sort(names(res))]

  return(res)
}


#' FORCE \eqn{K}-means solver.
#'
#' Solves a K-means SDP Relaxation using the FORCE algorithm when \eqn{K} is unknown.
#'
#' @param D a matrix \eqn{D} as defined above.
#' @param force_opts tuning parameters. \code{NULL} signifies defaults will be used.
#' @param D_Kmeans matrix to be used for initial integer solution. \code{NULL} signifies that \code{D} will be used.
#' @param X0 initial iterate. \code{NULL} signifies that it will be generated randomly from \code{D_Kmeans}. If supplied, \code{E} must be supplied as well.
#'
#' @return An object with following components
#' \describe{
#' \item{\code{Z_T}}{Final iterate of the projected gradient descent algorithm run on the smoothed eigenvalue problem.}
#' \item{\code{B_Z_T}}{Projection of \code{Z_T} to the border of the positive semi-definite cone.}
#' \item{\code{B_Z_T_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_Z_T}.}
#' \item{\code{Z_best}}{Iterate with best objective value found during projected gradient descent on the smoothed eigenvalue problem.}
#' \item{\code{B_Z_best}}{Projection of \code{Z_best} to the border of the positive semi-definite cone.}
#' \item{\code{B_Z_best_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_Z_T}.}
#' \item{\code{km_best}}{Best clustering in terms of objective value of the SDP relaxation. This is found by running a complete linkage clustering algorithm on the rows of the iterate \code{Z_t}.}
#' \item{\code{B_km}}{Partnership matrix corresponding to \code{km_best}.}
#' \item{\code{km_opt_val}}{Objective value of the \eqn{K}-means SDP relaxation at \code{B_km}.}
#' \item{\code{km_best_time}}{Time elapsed (in seconds) until \code{km_best} was found.}
#' \item{\code{dual_certified}}{1 if a dual certificate was found, and 0 otherwise.}
#' \item{\code{dual_certified_grad_iter}}{Number of gradient updates performed before a dual certificate was found.}
#' \item{\code{dual_certified_time}}{Time elapsed (in seconds) until dual certificate was found for \code{B_km}.}
#' \item{\code{grad_iter_best}}{Gradient iteration where \code{Z_best} was computed.}
#' \item{\code{grad_iter_best_time}}{Time elapsed (in seconds) when the update \code{grad_iter_best} was performed.}
#' \item{\code{total_time}}{Total time elapsed (in seconds) during call to \code{gforce.FORCE}.}
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
#' res <- gforce.FORCE_adapt(D)
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#' @references J. Peng and Y. Wei. Approximating K-means-type Clustering via Semidefinite Programming. \emph{SIAM Journal on Optimization}, 2007.
#' @references J. Renegar.  Efficient first-order methods for linear programming and semidefinite programming. \emph{arXiv:1409.5832}, 2014.
#'
#' @seealso \code{\link{gforce.defaults}}
#' @useDynLib GFORCE FORCE_adapt_R
#' @export
gforce.FORCE_adapt <- function(D,force_opts = NULL,D_Kmeans = NULL, X0 = NULL) {
  d <- ncol(D)
  # initialize E and ESI
  E <- 0.5*diag(d) + (1/(2*d))*matrix(1,ncol=d,nrow=d)
  E_EVEV <- eigen(E)
  E_V <- E_EVEV$vectors
  E_D <- diag(E_EVEV$values)
  E_sqrt <- E_V%*%(E_D^(0.5))%*%t(E_V)
  ESI <- solve(E_sqrt)

  if(is.null(force_opts)){
      force_opts <- gforce.defaults(d)
  }

  if(is.null(D_Kmeans)) {
      D_Kmeans <- D
  }

  if(is.null(X0)){
    X0 <- FORCE_adapt_init(D,E,force_opts$adapt_init_mode,force_opts$initial_mixing)
  }

  # make sure that initialization is valid
  if(sum(X0*D) >= sum(E*D)) {
    stop(sprintf('gforce.FORCE_adapt -- D^TX_0 >= D^TE. Check D, and if D is well specified, pass valid X0 as argument. %3.3f >= %3.3f',
          sum(X0*D),sum(E*D)))
  }

  res <- NULL

  # SHOULD HAVE CHECK TO MAKE SURE THAT OBJECTIVE VALUES OKAY
  C_result <- .C(FORCE_adapt_R,
          D = as.double(D),
          D_Kmeans = as.double(D_Kmeans),
          E = as.double(E),
          ESI = as.double(ESI),
          X0 = as.double(X0),
          d = as.integer(d),
          verbosity = as.integer(force_opts$verbose),
          kmeans_iter = as.integer(force_opts$kmeans_iter),
          dual_frequency = as.integer(force_opts$dual_frequency),
          max_iter = as.integer(force_opts$max_iter),
          finish_pgd = as.integer(force_opts$finish_pgd),
          primal_only = as.integer(force_opts$primal_only),
          number_restarts = as.integer(length(force_opts$restarts)),
          restarts = as.integer(force_opts$restarts),
          alpha = as.double(force_opts$alpha),
          eps_obj = as.double(force_opts$eps_obj),
          early_stop_mode = as.integer(force_opts$early_stop_mode),
          early_stop_lag = as.integer(force_opts$early_stop_lag),
          early_stop_eps = as.double(force_opts$early_stop_eps),
          Z_T = numeric(d^2),
          B_Z_T = numeric(d^2),
          Z_T_lmin = as.double(1.0),
          Z_best = numeric(d^2),
          B_Z_best = numeric(d^2),
          Z_best_lmin = as.double(2.0),
          B_Z_T_opt_val = as.double(3.0),
          B_Z_best_opt_val = as.double(4.0),
          km_opt_val = as.double(5.0),
          km_best = as.integer(numeric(d)),
          km_best_time = as.double(6.0),
          km_iter_best = as.integer(7),
          km_iter_total = as.integer(8),
          dc = as.integer(9),
          dc_time = as.double(10.0),
          dc_grad_iter = as.integer(11),
          grad_iter_best = as.integer(12),
          grad_iter_best_time = as.double(13.0),
          total_time = as.double(14))

  # Build Result
  res$Z_T <- matrix(C_result$Z_T,ncol=d)
  res$B_Z_T <- matrix(C_result$B_Z_T,ncol=d)
  res$B_Z_T_opt_val <- C_result$B_Z_T_opt_val
  res$Z_best <- matrix(C_result$Z_best,ncol=d)
  res$B_Z_best <- matrix(C_result$B_Z_best,ncol=d)
  res$B_Z_best_opt_val <- C_result$B_Z_best_opt_val
  res$km_best <- C_result$km_best
  res$B_km <- gforce.clust2mat(res$km_best)
  res$km_opt_val <- C_result$km_opt_val
  res$km_best_time <- C_result$km_best_time
  res$dual_certified <- C_result$dc
  res$dual_certified_grad_iter <- C_result$dc_grad_iter
  res$dual_certified_time <- C_result$dc_time
  res$grad_iter_best <- C_result$grad_iter_best
  res$grad_iter_best_time <- C_result$grad_iter_best_time
  res$total_time <- C_result$total_time

  res$X0 <- X0
  res$E <- E

  # perform HC clustering on final result
  # hc_res <- gforce.hclust(res$B_Z_best)
  # res$hc_best <- hc_res$clusters
  # res$hc_K <- hc_res$K
  res <- res[sort(names(res))]

  return(res)
}


#' Solve PECOK with FORCE.
#'
#' Uses the FORCE algorithm to solve the PECOK SDP.
#'
#' @param K number of clusters.
#' @param X \eqn{n x d} matrix. Either this or \code{D} must be specified.
#' @param D \eqn{d x d} matrix. Either this or \code{X} must be specified.
#' @param sigma_hat \eqn{d x d} matrix. If \code{D} is specified, this argument should be the
#' estimated covariance matrix. It is not strictly necessary to provide it, but it should be for
#' optimal performance. If \code{X} is specified, it will be ignored.
#' @param E strictly feasible solutions. \code{NULL} signifies that it will be generated randomly. If supplied, \code{X0} must be supplied as well.
#' @param gamma_par logical expression. If \code{gamma_par==TRUE}, then if \eqn{\Gamma} is computed,
#' a multi-threaded method is called, otherwise a single-threaded method is called.
#' @inheritParams gforce.FORCE
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#' @references J. Peng and Y. Wei. Approximating K-means-type Clustering via Semidefinite Programming. \emph{SIAM Journal on Optimization}, 2007.
#' @references F. Bunea, C. Giraud, M. Royer and N. Verzelen. PECOK: a convex optimization approach to variable clustering. \emph{arXiv:1606.05100}, 2016.
#'
#' @seealso \code{\link{gforce.defaults}}
#' @export
gforce.PECOK <- function(K, X=NULL, D=NULL, sigma_hat = NULL, force_opts = NULL, X0 = NULL, E = NULL, gamma_par = FALSE) {
  if(is.null(X) && is.null(D)) {
    stop('gforce.PECOK -- You must specify one of X or D.')
  } else if (!is.null(X) && !is.null(D)) {
    stop('gforce.PECOK -- You must specify one of X or D.')
  }
  if(is.null(D)){
    n <- nrow(X)
    gamma_hat <- gforce.Gamma(X)
    sigma_hat <- t(X)%*%X / n
    D <- diag(gamma_hat)-sigma_hat
  }
  if(is.null(sigma_hat)){
    sigma_hat <- D
  }

  res <- gforce.FORCE(D,K,D_Kmeans = sigma_hat, force_opts = force_opts, X0 = X0, E = E)
  res$D <- D
  return(res)
}



#' Solve PECOK Adaptive SDP with FORCE.
#'
#' Uses the FORCE algorithm to solve the PECOK SDP when \eqn{K} is unknown.
#'
#' @param X \eqn{n x d} matrix. Either this or \code{D} must be specified.
#' @param D \eqn{d x d} matrix. Either this or \code{X} must be specified.
#' @param sigma_hat \eqn{d x d} matrix. If \code{D} is specified, this argument should be the
#' estimated covariance matrix. It is not strictly necessary to provide it, but it should be for
#' optimal performance. If \code{X} is specified, it will be ignored.
#' @param gamma_par logical expression. If \code{gamma_par==TRUE}, then if \eqn{\Gamma} is computed,
#' a multi-threaded method is called, otherwise a single-threaded method is called.
#' @inheritParams gforce.FORCE_adapt
#'
#' @references C. Eisenach and H. Liu. Efficient, Certifiably Optimal High-Dimensional Clustering. \emph{arXiv:1806.00530}, 2018.
#' @references J. Peng and Y. Wei. Approximating K-means-type Clustering via Semidefinite Programming. \emph{SIAM Journal on Optimization}, 2007.
#' @references F. Bunea, C. Giraud, M. Royer and N. Verzelen. PECOK: a convex optimization approach to variable clustering. \emph{arXiv:1606.05100}, 2016.
#'
#' @seealso \code{\link{gforce.defaults}}
#' @export
gforce.PECOK_adapt <- function(X=NULL, D=NULL, sigma_hat = NULL, force_opts = NULL, X0 = NULL, gamma_par = FALSE) {
  if(is.null(X) && is.null(D)) {
    stop('gforce.PECOK -- You must specify one of X or D.')
  } else if (!is.null(X) && !is.null(D)) {
    stop('gforce.PECOK -- You must specify one of X or D.')
  }
  if(is.null(D)){
    n <- nrow(X)
    d <- ncol(X)
    gamma_hat <- gforce.Gamma(X)
    sigma_hat <- t(X)%*%X / n
    D_orig <- diag(gamma_hat)-sigma_hat
    kappa_hat <- max(gamma_hat)*(d/n + sqrt(d/n))
    D <- D_orig + kappa_hat*diag(d)
  }
  if(is.null(sigma_hat)){
    sigma_hat <- D
  }


  res <- gforce.FORCE_adapt(D,D_Kmeans = sigma_hat, force_opts = force_opts, X0 = X0)
  res$D <- D
  return(res)
}




#' FORCE default tuning parameters.
#'
#' Provides the default tuning parameters for \code{\link{gforce.FORCE}}.
#' @param d dimension of random vector or number of datapoints.
#'
#' @return An object with following components
#' \describe{
#' \item{\code{adapt_init_mode}}{a numeric. Indicates which initialization mode to use for \code{\link{gforce.FORCE_adapt}}.}
#' \item{\code{alpha}}{a numeric. Gives the step size for the projected gradient descent updates.}
#' \item{\code{dual_frequency}}{an integer. Specifies how many gradient updates to perform between searches for a dual certificate. }
#' \item{\code{duality_gap}}{a numeric. If the duality gap can be shown to be less than \code{duality_gap}, the FORCE algorithm terminates.}
#' \item{\code{early_stop_mode}}{a numeric. \code{early_stop_mode == 1} indicates that the algorithm should use an early stopping rule.}
#' \item{\code{early_stop_lag}}{an integer. This indicates the number of iterations without sufficient improvement in objective value before early stopping.}
#' \item{\code{early_stop_eps}}{a numeric. Threshold for objective value improvement used to determine early stopping.}
#' \item{\code{eps_obj}}{a numeric. Specifies the precision required of the optimal solution to the eigenvalue maximization problem.}
#' \item{\code{finish_pgd}}{an integer. If \code{finish_pgd} is 1, then other stopping criteria are ignored and FORCE performs \code{max_iter} gradient updates.}
#' \item{\code{initial_mixing}}{a numeric between 0 and 1. Specifies how to construct the initial strictly feasible solution to the SDP relaxation.}
#' \item{\code{kmeans_iter}}{an integer. The number of times to run a \eqn{K}-means solver during each search for an optimal clustering and dual certificate.}
#' \item{\code{max_iter}}{an integer. The maximum number of gradient updates to perform.}
#' \item{\code{primal_only}}{an integer. \code{primal_only == 1} indicates that the algorithm should not search for a dual certificate.}
#' \item{\code{restarts}}{a vector of integers. This specifies the iterations at which to take the projection of the current iterate and restart the algorithm with that as the initial solution.}
#' \item{\code{verbose}}{an integer. Specifies the level of verbosity requested from gforce.FORCE.}
#' }
#'
#' @examples
#' opts <- gforce.defaults(20)
#' @export
gforce.defaults <- function(d){
  options <- NULL
  options$adapt_init_mode = 0
  options$alpha = 10^-4
  options$dual_frequency = 50
  options$duality_gap = 10^-5
  options$early_stop_mode = 1
  options$early_stop_lag = 50
  options$early_stop_eps = 10^-5
  options$eps_obj = 0.01
  options$finish_pgd = 0
  options$initial_mixing = 2/d
  options$kmeans_iter = 10
  options$max_iter = 500
  options$primal_only = 0
  options$restarts = c(75)
  options$verbose = 0

  return(options)
}


FORCE_adapt_init <- function(D,E,mode,mixing) {
  X0 <- NULL
  d <- nrow(D)
  if(mode == 1) {
    hc_res <- gforce.hclust(D)
    hc_sol <- gforce.clust2mat(hc_res$clusters)
    X0 <- mixing*hc_sol + (1-mixing)*E
  } else if(mode==2) {
    hc_res <- gforce.hclust(D)
    hc_sol <- gforce.clust2mat(hc_res$clusters)
    X0 <- 0.5*hc_sol + 0.5*(matrix(1,ncol=d,nrow=d))*(1/d)
  } else if (mode==3) {
    X0 <- (1/d)*matrix(1,ncol=d,nrow=d)
  } else {
    hc_res <- gforce.hclust(D)
    hc_sol <- gforce.clust2mat(hc_res$clusters)
    X0_a <- mixing*hc_sol + (1-mixing)*E
    X0_b <- 0.5*hc_sol + 0.5*(matrix(1,ncol=d,nrow=d))*(1/d)
    val_a <- sum(X0_a*D)
    val_b <- sum(X0_b*D)
    X0 <- X0_a
    if(val_b < val_a) {
      X0 <- X0_b
    }
  }

  return(X0)
}
