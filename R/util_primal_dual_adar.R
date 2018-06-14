# primal dual method for the smoothed problem with adaptive restart
# see Candes, et al 2012

# solves MINIMIZATION form of SDP
#note - D is the D in the minimization form...not the same as D in
#PECOK SDP

primal_dual_smoothed_adar <- function(D,sigma,K,options,X0,E,ESI){
  if(options$verbose > -1){
    cat(sprintf("\tSolving K-Means SDP with FORCE\r\n"))
    cat(sprintf("\t\tOptions -- Verbosity: %d\r\n",options$verbose))
    cat(sprintf("\t\tOptions -- Finish PGD: %d\r\n",options$finish_pgd))
    cat(sprintf("\t\tOptions -- Max. Iter.: %d\r\n",options$max_iter))
    cat(sprintf("\t\tOptions -- Dual Frequency: %d\r\n",options$dual_frequency))
    cat(sprintf("\t\tOptions -- K-means Reps.: %d\r\n",options$kmeans_iter))
    cat(sprintf("\t\tOptions -- Eps. Objective: %.3f\r\n",options$eps_obj))
    cat(sprintf("\t\tOptions -- Alpha: %.5f\r\n",options$alpha))
  }

  start_time <- proc.time()
  # initialize return values structure fields to -1
  results <- NULL
  results$B_km <- -1  #set when found
  results$km_opt_val <- -1  #set when found
  results$B_Z_best <- -1  #set when_found
  results$dual_certified <- FALSE  #set when found
  results$dual_certified_grad_iter <- -1 
  results$dual_certified_time <- -1 
  results$dual_value_min <- -1 
  results$relative_error_bound <- 1 
  results$grad_iter_best <- -1  #set when found
  results$grad_iter_best_time <- -1 
  results$km_iter_best <- -1  #set when found
  results$km_best_time <- -1  #set when found
  grad_iter_total <- -1  #set as found
  results$km_iter_total <- -1  #set when found

  #initialize - get initial best clustering via kmeans
  results$km_best <- gforce.kmeans(sigma,K)$clusters
  results$B_km <- gforce.clust2mat(results$km_best)
  results$km_iter_best <- 1
  results$km_iter_total <- 1
  results$km_best_time <- proc.time() - start_time
  results$km_opt_val <- sum(-D*results$B_km)
  for(i in 1:(options$kmeans_iter-1)){
    results$km_iter_total <- results$km_iter_total + 1
    estimate_i_idx <- gforce.kmeans(sigma,K)$clusters

    estimate_i <- gforce.clust2mat(estimate_i_idx)
    estimate_i_val <- sum(-D*estimate_i)
    if(estimate_i_val > results$km_opt_val) {
      results$km_best <- estimate_i_idx
      results$B_km <- estimate_i
      results$km_iter_best <- results$km_iter_total
      results$km_opt_val <- estimate_i_val
      results$km_best_time <- proc.time() - start_time
    }
  }
  
  #initialize - PGD parameters
  p <- dim(D)[1]
  mu <- 0.5*options$eps_obj/log(p)
  alpha <- options$alpha
  
  #initialize - PGD initial values
  grad_iter_total <- 1
  results$grad_iter_best <- grad_iter_total
  results$grad_iter_best_time <- proc.time() - start_time
  outer_iterations <- 0
  X_tp1 <- X0 # == y_1 == z_1
  Z_tp1 <- X0 # == y_1 == z_1
  Z_best <- X0
  lambda_t <- 0 #initialized to lambda_t <- lambda_0
  lambda_tp1 <- 1 #initialized to lambda_{t+1} <- lambda_1
  s_res <- smoothed_objective(Z_tp1,E,ESI,mu)
  obj_tp1 <- s_res$objective_value
  lambda_min_tp1 <- s_res$lambda_min
  obj_best <- obj_tp1
  lambda_min_best <- lambda_min_tp1
  stop_criterion <- 1
  if(length(options$restarts) > 0){
    next_restart_t <- options$restarts[1]
    current_restart <- 1
  } else {
    next_restart_t <- -1
    current_restart <- -1
  }
  num_momentum_restarts <- 0

  if(options$verbose > 0){
    cat(sprintf("\tINITIAL VALUES - lambda_min_E(X_1) = %.4f\r\n",lambda_min_tp1))
    cat(sprintf("\tINITIAL VALUES - f_mu,E(X_1) = %.4f\r\n",obj_tp1))
  }

  #initialize - find initial dual solution
  kmeans_dual_feasible <- FALSE
  dual_s <- dual_solution(results$km_best,-D)
  if(!is.null(dual_s)){
    kmeans_dual_feasible <- TRUE
    results$dual_value_min <- dual_s$dual_value
    results$relative_error_bound <- (results$dual_value_min  - results$km_opt_val)/results$dual_value_min
    if(results$relative_error_bound < options$eps_obj){
      results$dual_certified <- TRUE
      results$dual_certified_grad_iter <- grad_iter_total
      results$dual_certified_time <- proc.time() - start_time
    }
  }
  


  #outer iterator
  while(stop_criterion && (results$relative_error_bound > options$eps_obj || !kmeans_dual_feasible || options$finish_pgd)){
    next_iter <- 1 #allows entry into inner iteration
    if(options$verbose > 0){
      cat(sprintf("\tOUTER ITERATION %d -- START\r\n",outer_iterations))
    }
    #inner iterator (gradient descent)
    while(((grad_iter_total %% options$dual_frequency) > 0 || next_iter) && stop_criterion){
      #check if need to restart the method
      if(options$verbose > 1){
        cat(sprintf("\t\tINNER ITERATION %d -- START\r\n",grad_iter_total))
      }

      if(grad_iter_total == next_restart_t){
        if(options$verbose > 1){
          cat(sprintf("\t\tINNER ITERATION %d -- RESTART X0\r\n",grad_iter_total))
        }
        X_tp1 <- E + ( 1/(1 - lambda_min_best) ) * (Z_best - E)
        Z_tp1 <- X_tp1
        s_res <- smoothed_objective(Z_tp1,E,ESI,mu)
        obj_tp1 <- s_res$objective_value
        lambda_min_tp1 <- s_res$lambda_min
        results$grad_iter_best <- grad_iter_total
        results$grad_iter_best_time <- proc.time() - start_time
        Z_best <- Z_tp1
        obj_best <- obj_tp1
        lambda_min_best <- lambda_min_tp1
        lambda_t <- 0
        lambda_tp1 <- 1
        if(current_restart < length(options$restarts)){
          current_restart <- current_restart + 1
          next_restart_t <- options$restarts[current_restart]
        }
      }
      
      #update - previous values
      next_iter <- 0
      X_t <- X_tp1
      Z_t <- Z_tp1
      obj_t <- obj_tp1
      lambda_min_t <- lambda_min_tp1
      
      #update - auxiliary sequences
      lambda_t <- lambda_tp1
      lambda_tp1 <- (1+ (1 + 4*lambda_t^2)^0.5)/2
      gamma_t <- (1 - lambda_t) / lambda_tp1
      
      #update - find gradient
      grad_t <- smoothed_gradient(X_t,E,ESI,mu)
      GX_t <- grad_t$GX
      GS_t <- grad_t$GS

      #update - primary sequences
      c_perp_res <- C_perp_update(X_t,Z_t,obj_t,GX_t,GS_t,E,ESI,mu,D,alpha,options)
      Z_tp1 <- c_perp_res$Z_tp1
      obj_tp1 <- c_perp_res$obj_tp1
      lambda_min_tp1 <- c_perp_res$lambda_min_tp1
      X_tp1 <- (1 - gamma_t)*Z_tp1 + gamma_t*Z_t
      
      #output status

      if(options$verbose > 1) {
        proj_Z_tp1 <- E + ( 1/(1 - lambda_min_tp1) ) * (Z_tp1 - E)
        cat(sprintf('\t\tINNER ITERATION %d -- SDP Objective Value: %f\r\n',grad_iter_total,sum(-D*proj_Z_tp1)))
        cat(sprintf('\t\tINNER ITERATION %d -- f_mu(Z_{t+1}) <- %f\r\n',grad_iter_total,obj_tp1))
        cat(sprintf('\t\tINNER ITERATION %d -- lambda_min_E(Z_{t+1}) <- %f\r\n',grad_iter_total,lambda_min_tp1))
        cat(sprintf('\t\tINNER ITERATION %d -- max_{ij} |(GX_t)_{ij}| <- %f\r\n',grad_iter_total,max(max(abs(GX_t)))))
        cat(sprintf('\t\tINNER ITERATION %d -- max_{ij} |(GS_t)_{ij}| <- %f\r\n',grad_iter_total,max(max(abs(GS_t)))))
        cat(sprintf('\t\tINNER ITERATION %d -- ||GX_t||_F <- %f\r\n',grad_iter_total,norm(GX_t,'F')))
        cat(sprintf('\t\tINNER ITERATION %d -- ||GS_t||_F <- %f\r\n',grad_iter_total,norm(GS_t,'F')))
      }
      #check to see if need to restart momentum term
      if(obj_tp1 < obj_t){
        # restart
        lambda_t <- 0 #initialized to lambda_t <- lambda_0
        lambda_tp1 <- 1 #initialized to lambda_{t+1} <- lambda_1
        if(options$verbose > 1) {
          cat(sprintf('\t\tINNER ITERATION %d -- Restarting Momentum\r\n',grad_iter_total))
        }
        num_momentum_restarts <- num_momentum_restarts + 1
      }

      #update counter, best candidate solution found
      if(obj_tp1 > obj_best){
        if(options$verbose > 1) {
          cat(sprintf('\t\tINNER ITERATION %d -- Found New Best\r\n',grad_iter_total))
        }
        obj_best <- obj_tp1
        Z_best <- Z_tp1
        lambda_min_best <- lambda_min_tp1
        results$grad_iter_best <- grad_iter_total + 1
        results$grad_iter_best_time <- proc.time() - start_time
      }
      if(options$verbose > 1) {
        cat(sprintf('\t\tINNER ITERATION %d -- COMPLETE\r\n',grad_iter_total))
      }
      grad_iter_total <- grad_iter_total + 1
      # check the stopping criterion
      stop_criterion <- check_stopping_criterion(grad_iter_total,results$grad_iter_best,options)            
    }
    
    results$total_time <- proc.time() - start_time

    #find better (hopefully) candidate integer solution
    proj_Z_tp1 <- E + ( 1/(1 - lambda_min_tp1) ) * (Z_tp1 - E)
    found_better_kmeans <- 0
    for(i in 1:(options$kmeans_iter)){
      estimate_i_idx <- gforce.kmeans(proj_Z_tp1,K)$clusters
      estimate_i <- gforce.clust2mat(estimate_i_idx)
      estimate_i_val <- sum(-D*estimate_i)
      results$km_iter_total <- results$km_iter_total + 1
      if(estimate_i_val > results$km_opt_val){
        results$km_iter_best <- results$km_iter_total
        found_better_kmeans <- 1
        results$km_best <- estimate_i_idx
        results$B_km <- estimate_i
        results$km_opt_val <- estimate_i_val
        results$km_best_time <- proc.time() - start_time
      }
    }
    if(found_better_kmeans) {
      dual_s <- dual_solution(results$km_best,-D)
      if(!is.null(dual_s)){
        results$dual_value_min <- dual_s$dual_value
        results$relative_error_bound <- (results$dual_value_min  - results$km_opt_val)/results$dual_value_min
        kmeans_dual_feasible <- TRUE
      }
      
    }
    if(found_better_kmeans && results$relative_error_bound < options$eps_obj){
      results$dual_certified <- 1
      results$dual_certified_grad_iter <- grad_iter_total
      results$dual_certified_time <- proc.time() - start_time
    }
    outer_iterations <- outer_iterations + 1
    if(options$verbose > 0){
      cat(sprintf('\t\tOUTER ITERATION %d -- Dual Certified: %d\r\n',outer_iterations,results$dual_certified))
      cat(sprintf('\t\tOUTER ITERATION %d -- COMPLETE\r\n',outer_iterations))
    }
  }

  results$Z_T <- Z_tp1
  results$B_Z_T <- E + ( 1/(1 - lambda_min_tp1) ) * (Z_tp1 - E)
  results$B_Z_T_opt_val <- sum(-D*results$B_Z_T)
  results$Z_best <- Z_best
  results$B_Z_best <- E + ( 1/(1 - lambda_min_best) ) * (Z_best - E)
  results$B_Z_best_opt_val <- sum(-D*results$B_Z_best)

  #print results
  if(options$verbose > -1){
      cat(sprintf("\tFORCE Algorithm Complete\r\n"))
      cat(sprintf("\t\tDual Certificate:%d\r\n",results$dual_certified))
      cat(sprintf("\t\tInner Iterations Total:%d\r\n",results$grad_iter_total))
      cat(sprintf("\t\tOuter Iterations Total:%d\r\n",outer_iterations))
      cat(sprintf("\t\tTotal Running Time:%.3fs\r\n",0))
      cat(sprintf("\t\t<D,B_Z_T> = %.4f\r\n",results$B_Z_T_opt_val))
      cat(sprintf("\t\t<D,B_Z_best> = %.4f\r\n",results$B_Z_best_opt_val))
      cat(sprintf("\t\t<D,B_km> = %.4f\r\n",results$km_opt_val))
  }
  if(options$verbose == 5){
      cat(sprintf("\t\tTotal Momentum Restarts:%d\r\n",num_momentum_restarts))
  }

  return(results)
}

check_stopping_criterion <- function(t,t_best,options){
  contin <- TRUE
  if(t > options$max_iter){
    contin <- FALSE
  }
  return(contin)
}

C_perp_update <- function(X_t,Z_t,obj_t,GX_t,GS_t,E,ESI,mu,D,alpha,options){
  grad_perp <- project_C_perpendicular(GX_t,GS_t,D,options$slack_scale)
  G_t_perp <- grad_perp$Z_proj
  norm_G_t_perp <- 2*sqrt(1 + 1/options$slack_scale)*norm(G_t_perp,'F') # need to mult by two for slacks
  res <- NULL

  G_t_perp <- G_t_perp / norm_G_t_perp        
  Z_tp1 <- X_t + alpha*G_t_perp
  s_res <- smoothed_objective(Z_tp1,E,ESI,mu)
  res$obj_tp1 <- s_res$objective_value
  res$lambda_min_tp1 <- s_res$lambda_min
  res$Z_tp1 <- Z_tp1

  # if(options$verbose > 0){
  #   corr <- 2*sum(G_t_perp*(Z_t - X_t))
  #   alpha_criterion <- 2*(res$obj_tp1 - obj_t + corr)/(norm_G_t_perp^2)
  #   cat(sprintf('\t\tProjection -- Norm of G_t perpendicular -- %f\r\n',norm_G_t_perp))
  #   cat(sprintf('\t\tProjection -- Correlation of gradient, projection -- %f\r\n',corr))
  #   cat(sprintf('\t\tProjection -- Alpha Criterion -- %f\r\n',alpha_criterion))
  #   cat(sprintf('\t\tProjection -- alpha -- %f\r\n',alpha))
  # }
  return(res)
}
