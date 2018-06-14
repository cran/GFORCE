
#objective function takes only two argunments -- data and parameter
#solver only takes two arguments -- data and tuning parameter
cv_lambda_selection <- function(objective_function,function_solver,X_vals,group_assignments,
                                  lambda_min_exp,lambda_max_exp,num_folds,subdivisions=3,graph='latent'){
  # create test and train data
  n <- nrow(X_vals)
  groups <- unique(group_assignments)
  K <- length(groups)
  fold_assignments <- cross_validation_folds(n,num_folds)
  
  #compute C hat or S_hat on all folds 
  train_data <- array(0,c(K,K,num_folds))
  test_data <- array(0,c(K,K,num_folds))
  if(graph=='latent') {
    for(i in 1:num_folds){
      X <- X_vals[fold_assignments == i,] 
      sig_hat <- (t(X)%*%X)/n
      gam_hat <- glatent_Gamma_hat(sig_hat,group_assignments)
      test_data[,,i] <- C_hat(sig_hat,gam_hat,group_assignments)
      X <- X_vals[fold_assignments != i,]
      sig_hat <- (t(X)%*%X)/n
      gam_hat <- glatent_Gamma_hat(sig_hat,group_assignments)
      train_data[,,i] <- C_hat(sig_hat,gam_hat,group_assignments)
    }
  } else {
    for(i in 1:num_folds){
      X <- X_vals[fold_assignments == i,] 
      sig_hat <- (t(X)%*%X)/n
      test_data[,,i] <- S_hat(sig_hat,group_assignments)
      X <- X_vals[fold_assignments != i,]
      sig_hat <- (t(X)%*%X)/n
      train_data[,,i] <- S_hat(sig_hat,group_assignments)
    }
  }

  # setup -- get function pointer and generate lambda grid
  lambda_sequence <- 2^seq(lambda_min_exp,lambda_max_exp,by=1/subdivisions)
  lambda_total <- length(lambda_sequence)
  lambda_test_values <- rep(0,lambda_total)
  num_folds <- dim(train_data)[3]
  
  # try each lambda in grid
  for(i in 1:lambda_total){
    lambda <- lambda_sequence[i] 
    for(j in 1:num_folds){
      train_estimator <- function_solver(train_data[,,j],lambda) 
      lambda_test_values[i] <- lambda_test_values[i] + objective_function(test_data[,,j],train_estimator)
    }
  }
  return(lambda_sequence[which.min(lambda_test_values)])
}

# returns indices for num-folds-fold cross-validation sample split
cross_validation_folds <- function(num_samples,num_folds){
  shuffled_idx <- sample(num_samples)
  fold_size_smaller <- floor(num_samples/num_folds)
  fold_size_larger <- fold_size_smaller + 1
  num_larger_folds <- num_samples %% num_folds
  fold_idx <- rep(0,num_samples)
  cur_idx <- 1
  if(num_larger_folds > 0){
    for(i in 1:num_larger_folds){
      last_idx <- cur_idx + fold_size_larger - 1
      fold_idx[shuffled_idx[cur_idx:last_idx]] <- i
      cur_idx <- last_idx + 1
    }
  }
  cur_fold <- num_larger_folds + 1
  while(cur_fold <= num_folds){
    last_idx <- cur_idx + fold_size_smaller -1
    fold_idx[shuffled_idx[cur_idx:last_idx]] <- cur_fold
    cur_idx <- last_idx + 1
    cur_fold <- cur_fold + 1
  }
  return(fold_idx)
}