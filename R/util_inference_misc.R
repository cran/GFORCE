# Misc Utility Functions for GBlock Inference

coverage_percentage <- function(theta_star,conf_intervals){
  K <- nrow(theta_star)
  covered <- (theta_star > conf_intervals[,,1]) & (theta_star < conf_intervals[,,3])
  return(sum(colSums(covered))/K^2)
}


error_rates <- function(matrix_confs,pop_matrix){
  true_null <- 0
  false_positives <- 0
  false_negatives <- 0
  nan_found <- 0
  K <- ncol(pop_matrix)
  
  for(t in 1:K){
    for(k in 1:K){
      if(!is.na(matrix_confs[t,k,1])) {
        if(abs(pop_matrix[t,k]) > 10^(-6)){
          if(matrix_confs[t,k,1] <= 0 && matrix_confs[t,k,3] >= 0 ){
            #don't reject -- missed positive
            false_negatives <- false_negatives + 1
          }else{
            #correctly rejected
            #do nothing
          }
        } else{
          #true null
          true_null <- true_null + 1
          if(matrix_confs[t,k,1] <= 0 && matrix_confs[t,k,3] >= 0 ){
            #correctly fail to reject
          }else{
            #rejected incorrectly -- false positive
            false_positives <- false_positives + 1
          }
        }
      } else {
        nan_found <- nan_found + 1
      }
    }
  }
  
  true_positive <- K^2 - true_null - nan_found
  res <- NULL
  if(true_null == 0){
    res$type_1_rate <- 0
  } else {
    res$type_1_rate <- false_positives / true_null
  }
  if(true_positive == 0){
    res$type_2_rate <- 0
  } else {
    res$type_2_rate <- false_negatives / true_positive
  }
  
  lower_bound <- matrix_confs[,,1]
  upper_bound <- matrix_confs[,,3]
  lower_bound[is.na(lower_bound)] <- Inf
  upper_bound[is.na(upper_bound)] <- -Inf
  res$coverage_percentage <- sum((lower_bound <= pop_matrix) & (pop_matrix <= upper_bound))/(K^2)
  res$nan_found <- nan_found
  return(res)
}

order_group_assignments <- function(group_assignments){
  ordering <- unique(group_assignments)
  K <- length(ordering)
  new_group_assignments <- rep(0,length(group_assignments))
  for(i in 1:K){
    group_idx <- which(group_assignments == ordering[i])
    new_group_assignments[group_idx] <- i
  }
  return(new_group_assignments)
}

Xi_star <- function (cstar,group_assignments){
  group_assignments <- order_group_assignments(group_assignments)
  sig_star <- A_hat(group_assignments)%*%cstar%*%t(A_hat(group_assignments))
  s_star <- S_hat(sig_star,group_assignments)
  xi_star <- solve(s_star)
  return(xi_star)
}


# S hat estimator
S_hat <- function(sig_hat,group_assignments){
  Ah <- A_hat(group_assignments)
  ATA <- t(Ah)%*%Ah
  ATAI <- solve(ATA)
  Sh <- ATAI%*%t(Ah)%*%(sig_hat)%*%Ah%*%ATAI
  return(Sh)
}

C_hat <- function(sig_hat,gam_hat,group_assignments){
  Ah <- A_hat(group_assignments)
  ATA <- t(Ah)%*%Ah
  ATAI <- solve(ATA)
  Ch <- ATAI%*%t(Ah)%*%(sig_hat - diag(gam_hat))%*%Ah%*%ATAI
  return(Ch)
}

A_hat <- function(group_assignments) {
  d <- length(group_assignments)
  K <- length(unique(group_assignments))
  Ah <- matrix(rep(0,K*d),nrow=d)
  for(i in 1:d) {
    Ah[i,group_assignments[i]] <- 1
  }
  return(Ah)
}


# Implements \Hat\bGamma estimator
# WARNING  --- ASSUMES GROUPS NUMBERED 1..K
# ALSO Chat will be according to GROUP NAMES
# NOT PERMUTATION INVARIANT

#assumes group assignments 
glatent_Gamma_hat <- function(sig_hat,group_assignments) {
  gh <- vector(length = length(group_assignments))
  groups <- unique(group_assignments)
  K <- length(groups)
  
  for(k in 1:K) {
    group_idx <- which(group_assignments == k)
    group_size <- length(group_idx)
    if(group_size == 1){
      gh[group_idx] <- 0
    } else {
      sig_group <- sig_hat[group_idx,group_idx]
      sig_group_colsums <- colSums(sig_group)
      sig_group_vars <- diag(sig_group)
      gh[group_idx] <- (group_size/(group_size - 1))*sig_group_vars - (1/(group_size-1))*sig_group_colsums
    }
  }
  
  return(gh)
}

# Estimate \hat w
w_hat <- function(Chat,t,lambda) {
  Km1 <- nrow(Chat) - 1
  A_orig <- -1*Chat[-t,-t]
  #  print(A_orig)
  lp_rhs <- vector(length=4*Km1)
  
  A <- matrix(rep(0,8*Km1^2),nrow=4*Km1)
  rel <- rep("<=",4*Km1)
  obj <- c(rep(1,Km1),rep(0,Km1)) 
  
  # RHS of l_infty norm constraints
  lp_rhs[1:(2*Km1)] <- lambda*rep(1,2*Km1) + c(-1*Chat[t,-t],Chat[t,-t])

    # RHS of l_infty norm constraints
  lp_rhs[(2*Km1+1):(4*Km1)] <- rep(0,2*Km1)

  # l_infty norm constraints
  A[(1:Km1),1:Km1] <- A_orig
  A[(1:Km1),(Km1+1):(2*Km1)] <- -1*A_orig
  A[(Km1+1):(2*Km1),1:Km1] <- -1*A_orig
  A[(Km1+1):(2*Km1),(Km1+1):(2*Km1)] <- A_orig
  
  # solve
  res <- lpSolve::lp(direction="min",obj,A,rel,lp_rhs)

  return(res$solution[1:Km1] - res$solution[(Km1+1):(2*Km1)])
}
 
v_hat <- function(Chat,t,lambda) {
  wh <- w_hat(Chat,t,lambda)
  vt <- rep(1,length(wh)+1)
  vt[-t] = -1*wh
  return(vt)
}


latent_spectrum_check <- function(Chat,epsilon=0.00001) {
  edecomp <- eigen(Chat)
  evals <- edecomp$values
  evecs <- edecomp$vectors
  eval_shift_idx <- which(evals < epsilon)
  evals[eval_shift_idx] <- epsilon

  Chat_shift <- evecs%*%diag(evals)%*%t(evecs)

  return(Chat_shift)

  # edecomp$vectors%*%diag(edecomp$values)%*%t(edecomp$vectors)
}