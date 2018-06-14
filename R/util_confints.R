
#decorelated estimator
decorrelated_estimator <- function(C,theta_k,v,t,k){
  return( as.numeric((v[k] - v%*%C[,-t]%*%theta_k[-t]) / (v%*%C[,t]) ))
}

decorrelated_estimator_t <- function(C,theta_hat,v_t,t) {
  prod <- as.numeric(v_t%*%C[,-t]%*%theta_hat[-t,])
  return((v_t - prod) / as.vector(v_t%*%C[,t]))
}


#comnputes \bar{\Gamma} from Appendix D of GBlock Paper
variance_latent_group_average_error <- function(gh,group_assignments){
  K <- length(unique(group_assignments))
  gh_bar <- rep(0,K)
  for(i in 1:K){
    group_idx <- which(group_assignments == i)
    group_size <- length(group_idx)
    gh_bar[i] <- sum(gh[group_idx]) / (group_size^2)
  }
  
  return(gh_bar)  
}

# computes the diagonal matrix P in Appendix D of GBlock Paper
# GBlock assumes |G_i| > 1, but if G_i = 1 we can assume gamma_i = 0
# and the inference results still hold. In this case, the LOT in
# the variance estimate are just 0.
variance_latent_lot_part_B <- function(gh,group_assignments) {
  K <- length(unique(group_assignments))
  lot_mat <- rep(0,K)
  
  for(i in 1:K){
    group_idx <- which(group_assignments == i)
    group_size <- length(group_idx)
    cross_sums <- 0
    #check if group is singleton
    if(group_size <= 1){
      lot_mat[i] <- 0
    } else {
      for(l in 1:group_size){
        for(m in 1:group_size){
          if(m != l){
            cross_sums <- cross_sums + gh[group_idx[l]]*gh[group_idx[m]]
          }
        }
      }
      sq_sums <- sum(gh[group_idx]^2)
      gh_bar_i_sq <- sum(gh[group_idx])^2
      lot_mat[i] <- 8*sq_sums/(group_size^4) + 2*cross_sums/(group_size^2 * (group_size-1)^2) - gh_bar_i_sq/(group_size^4)
    }
  }
  
  return(lot_mat)
}

latent_graph_variances <- function(theta_hat,v_hat,gam_hat,group_assignments) {
  K <- ncol(theta_hat)

  latent_variances <- matrix(0,ncol=K,nrow=K)

  gh_bar_diag <- variance_latent_group_average_error(gam_hat,group_assignments)
  gh_bar_theta <- diag(gh_bar_diag)%*%theta_hat
  theta_gh_bar_theta <- t(theta_hat)%*%gh_bar_theta #only want diag, need to transpose to match up
  theta_gh_bar_theta <- diag(theta_gh_bar_theta)
  vt_gh_bar_theta <- v_hat%*%gh_bar_theta #need to index it t,k for correctness
  vt_gh_bar_vt <- v_hat%*%diag(gh_bar_diag)%*%t(v_hat)
  vt_gh_bar_vt <- diag(vt_gh_bar_vt)
  lot_part_B_diag <- variance_latent_lot_part_B(gam_hat,group_assignments)

  for(k in 1:K){
    lot_part_B_theta_k <- (theta_hat[,k]^2)*lot_part_B_diag

    vt_lot_part_B_theta_k_vt <- v_hat%*%diag(lot_part_B_theta_k)%*%t(v_hat)
    vt_lot_part_B_theta_k_vt <- diag(vt_lot_part_B_theta_k_vt)
    for(t in 1:K){
      variance_tk <- (theta_hat[t,k]^2 + theta_hat[t,t]*theta_hat[k,k])/(theta_hat[t,t]^2)
      variance_lot_part_A <- 2*theta_hat[t,k]*vt_gh_bar_theta[t,k]/theta_hat[t,t] + theta_gh_bar_theta[k]/theta_hat[t,t]
      variance_lot_part_A <- variance_lot_part_A + theta_hat[k,k]*vt_gh_bar_vt[t] + theta_gh_bar_theta[k]*vt_gh_bar_vt[t]
      variance_lot_part_B <- vt_lot_part_B_theta_k_vt[t]
      
      variance_tk <- (variance_tk + variance_lot_part_A+variance_lot_part_B)*(theta_hat[t,t]^2)
      latent_variances[t,k] <- variance_tk
    }
  }
  return(latent_variances)
}
