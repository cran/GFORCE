#' FDR Control Procedure.
#'
#' Performs by default the Benjamini-Yekutieli FDR control procedure. Optionally, the Benjamini-Hochberg thresholding
#' rule can be used instead. As input it takes a symmetric matrix of test statistics with standard normal null
#' distributions. The number of hypotheses tested is \eqn{d(d-1)/2}.
#' 
#' @param test_stats \eqn{d x d} symmetric matrix of test statistics.
#' @param alpha alpha level for the FDR control procedure.
#' @param procedure a string. \code{procedure == 'BY'} indicates to use the Benjamini-Yekutieli thresholding rule.
#'        \code{procedure == 'BH'} indicates to use the Benjamini-Hochberg thresholding rule.
#'
#' @return An object with following components
#' \describe{
#' \item{\code{reject_null}}{\eqn{d x d} upper triangular matrix. \code{TRUE} entries indicate the null hypothesis should be rejected.}
#' \item{\code{R_tau_hat}}{an integer. Indicates the number of hypotheses rejected.}
#' \item{\code{tau_hat}}{a real number. Indicates a threshold above which the null hypothesis is rejected.}
#' \item{\code{num_hypotheses}}{an integer. Indicates the number of hypotheses tested.}
#' }
#'
#' @references Y. Benjamini and Y. Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
#'             \emph{Journal of the Royal Statistical Society}, 1995.
#' @references Y. Benjamini and D. Yekutieli. The Control of the False Discovery Rate in Multiple Testing Under Dependency.
#'             \emph{The Annals of Statistics}, 2001.
#' @export
gforce.FDR_control <- function(test_stats,alpha,procedure='BY') {
    # set up
    K <- nrow(test_stats)
    num_hypotheses <- (K^2 - K)/2
    beta <- NULL
    if(procedure == 'BY') {
        N_BY <- sum(1/(1:num_hypotheses))
        beta <- alpha / (2*num_hypotheses*N_BY)
    } else if(procedure == 'BH') {
        beta <- alpha / (2*num_hypotheses)
    } else {
        stop('gforce.FDR_control -- procedure must be BH or BY.')
    }
    

    # create absolute values of test_stats for
    # procedure
    test_stats_step_up <- abs(test_stats)
    for(i in 1:K){
        for(j in 1:i) {
            test_stats_step_up[i,j] <- 0
        }
    }
    # get descending tau value possibilities
    tau_levels <- sort(test_stats_step_up,decreasing=TRUE)
    tau_levels <- tau_levels[1:num_hypotheses]

    # reject greater than or equal to tau
    R_tau_hat <- 0
    tau_hat <- Inf
    q_tau_hat <- Inf
    for(R_tau_new in 1:num_hypotheses){
        min_tau_new <- tau_levels[R_tau_new]
        p_tau <- beta*R_tau_new
        q_tau <- abs(stats::qnorm(p_tau))
        if(min_tau_new >= q_tau){
            tau_hat <- min_tau_new
            R_tau_hat <- R_tau_new
            q_tau_hat <- q_tau
        }
    }

    # true_discoveries <- test_stats_step_up > tau_hat
    true_discoveries <- test_stats_step_up >= tau_hat

    res <- NULL
    res$reject_null <- true_discoveries
    res$R_tau_hat <- R_tau_hat
    res$tau_hat <- tau_hat
    res$num_hypotheses <- num_hypotheses

    return(res)
}


#' Convert confidence intervals to equivalent test statistics.
#'
#' Can convert a 4D array encoding the confidence intervals for a precision
#' matrix to standard normal test-statistics.
#' 
#' @param conf_ints \eqn{d x d x 3} array. Each \eqn{d x d x 1} slice is a symmetric matrix.
#' @param alpha confidence level level of the confidence intervals.
#' 
#' @return a \eqn{d x d} symmetric matrix of test statistics. 
#'
#'
#' @export
gforce.confint2test <- function(conf_ints,alpha) {
    K <- nrow(conf_ints)
    test_stats <- matrix(rep(0,K^2),nrow=K)

    z_alpha <- stats::qnorm(1 - (alpha/2))

    for(i in 1:K){
        for(j in 1:K){
            v_ij <- conf_ints[i,j,3] - conf_ints[i,j,1]
            v_ij <- v_ij / (2*z_alpha)
            test_stats[i,j] <- conf_ints[i,j,2]/v_ij #/(conf_ints[i,i,2]*v_ij)
        }
    }

    return(test_stats)
}