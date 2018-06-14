
#' Estimates \eqn{\Gamma} for the PECOK SDP.
#' 
#' In particular, it returns in diagonal form the estimator
#' \eqn{\Gamma} used to construct the PECOK penalized covariance estimator.
#' @param X \eqn{n x d} matrix. Each row represents a realization of a
#' \eqn{d} dimensional random vector.
#' @param par logical expression. If \code{par == TRUE}, then a multi-threaded version
#' of the function is called. If \code{par == FALSE}, a single-threaded version is called.
#' @param fast_estimator logical expression. If \code{fast_estimator == TRUE}, then the alternative
#' estimator for \eqn{\hat \Gamma} is used.
#' @param R_only logical expression. If \code{R_only == TRUE}, then no native code is run. If 
#' \code{fast_estimator != TRUE} this is ignored.
#' @return The estimator \eqn{\Gamma} as a \eqn{d} dimensional numeric vector.
#'
#' @examples
#' dat <- gforce.generator(5,20,20,3)
#' gam_hat <- gforce.Gamma(dat$X)
#'
#' @references F. Bunea, C. Giraud, M. Royer and N. Verzelen. PECOK: a convex optimization approach to variable clustering. \emph{arXiv:1606.05100}, 2016.
#'
#' @useDynLib GFORCE v_measure
#' @useDynLib GFORCE v_measure_par
#' @useDynLib GFORCE gamma_alternative_estimator_R
#' @useDynLib GFORCE gamma_alternative_estimator_par_R
#' @export
gforce.Gamma <- function(X,par = FALSE,fast_estimator = FALSE,R_only=FALSE){
  dims <- dim(X)
  n <- dims[1]
  d <- dims[2]
  
  gd <- rep(0,d)
  ips <- t(X)%*%X
  ips_diag <- matrix(diag(ips),ncol=1)

  # use original or alternative estimator from PECOK paper
  if(!fast_estimator){
    ones_d <- matrix(rep(1,d),ncol=1)
    n_xc_xd <- (ones_d%*%t(ips_diag) + ips_diag%*%t(ones_d) - 2*ips)^0.5
    
    if(par){
      result <- .C("v_measure_par",
                   IPS=as.double(ips),
                   n_xc_xd=as.double(n_xc_xd),
                   dimension=as.integer(d),
                   vm=numeric(d^2))
    } else{
      result <- .C("v_measure",
                   IPS=as.double(ips),
                   n_xc_xd=as.double(n_xc_xd),
                   dimension=as.integer(d),
                   vm=numeric(d^2))
    }
    
    Vs_upper <- matrix(result$vm,ncol=d)
    Vs <- Vs_upper + t(Vs_upper)
    for(a in 1:d){
      # get neighbors
      v_a <- Vs[a,]
      vamax <- max(v_a)
      v_a[a] <- 2*vamax
      ne1 <- which.min(v_a)
      v_a[ne1] <- 2*vamax
      ne2 <- which.min(v_a)
      gd[a] <- (1/n)*(ips[a,a] + ips[ne1,ne2] - ips[a,ne1] - ips[a,ne2])
    }
  } else {
    ips_diag_sqrt <- ips_diag^0.5
    if(!R_only){
      if(par){
        result <- .C("gamma_alternative_estimator_par_R",
                     IPS=as.double(ips),
                     ips_diag_sqrt=as.double(ips_diag_sqrt),
                     dimension=as.integer(d),
                     scaling=as.double(1/n),
                     nes=as.integer(rep(0,d)),
                     gd=numeric(d))
      } else{
        result <- .C("gamma_alternative_estimator_R",
                     IPS=as.double(ips),
                     ips_diag_sqrt=as.double(ips_diag_sqrt),
                     dimension=as.integer(d),
                     scaling=as.double(1/n),
                     nes=as.integer(rep(0,d)),
                     gd=numeric(d))
      }
      gd <- result$gd
    } else {
      nes <- gamma_alternative_estimator_ne(ips,ips_diag_sqrt,d)
      # construct estimator
      for(a in 1:d){
        gd[a] <- (1/n)*(ips[a,a] - ips[a,nes[a]])
      }
    }
  }
  
  return(gd)
}

# R implementation of neighbor finding for alternative
# gamma estimator
gamma_alternative_estimator_ne <- function(ips,ips_diag_sqrt,d){
  nes <- rep(0,d)
  for(a in 1:d){
    min_val <- 0
    min_idx <- -1
    for(b in 1:d){
      if(a != b) {
        max_val <- -1
        for(c in 1:d){
          if(a != c && b != c){
            new_val <- (ips[a,c]-ips[b,c])/ips_diag_sqrt[c]
            new_val <- abs(new_val)
            if(new_val > max_val){
              max_val <- new_val
            }
          }
        }
        if(max_val < min_val || min_idx == -1){
          min_val <- max_val
          min_idx <- b
        }
      }
    }
    nes[a] <- min_idx
  }
  return(nes)
}