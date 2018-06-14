#include "math.h"
#include "stdio.h"
#include "time.h"
#include "R.h"
#include "R_ext/Lapack.h"
#include "FORCE.h"
#include "util_FORCE.h"
#include "util.h"
#include "util_mops.h"

// Needs for all pointers in return values to point to sufficiently allocated
// memory. Solves primal problem in MINIMIZATION FORM
// K IS UNKNOWN
// verbosity: -1 == all output suppressed
//             0 == startup info and termination info only
//             1 == startup, termination, iteration counter
//             2 == startup, termination, iteration full information
//             5 == DEBUG. Prints all info above plus some additional. It WILL trigger extra computation
void FORCE_adapt(double* D, double* D_kmeans, double* E, double* ESI,
                    double* X0, int d, pgd_opts* opts, pgd_results* results) {
    /////////////////////////////////////////////////////////////
    //// STEP 1 - Local Variable Initialization, Memory Allocation
    /////////////////////////////////////////////////////////////
    double dtmp1,dtmp2;
    double* ptmp1 = NULL;
    int d2 = d*d;
    double eps_obj = opts -> eps_obj;

    // Create Workspace, Initialize Problem Instance, Allocate Memory Pool
    problem_instance prob;
    workspace work;
    mem_pool free_d2;
    double mu = 0.5*eps_obj/log(d);
    initialize_problem_instance(D, E, ESI, mu, d, 0, &prob);
    allocate_workspace_FORCE_adapt(d, &work);
    free_d2.base = (void **) R_alloc(5,sizeof(void*));
    free_d2.length = 5;
    free_d2.start_idx=0;
    free_d2.end_idx=-1; //because pool starts empty
    ptmp1 = (void *) R_alloc(d2,sizeof(double));
    mem_pool_insert(&free_d2, ptmp1);
    ptmp1 = (void *) R_alloc(d2,sizeof(double));
    mem_pool_insert(&free_d2, ptmp1);
    ptmp1 = (void *) R_alloc(d2,sizeof(double));
    mem_pool_insert(&free_d2, ptmp1);
    ptmp1 = (void *) R_alloc(d2,sizeof(double));
    mem_pool_insert(&free_d2, ptmp1);

    // Non-convex rounding -- uses HC not KM
    hclust_t tmp_hc_sol;
    int* km_clusters_tmp = NULL;
    int* km_clusters_new = (int *) R_alloc(d,sizeof(int));
    int* km_clusters_best = (int *) R_alloc(d,sizeof(int));
    double km_val_best;
    double km_val_new;
    double km_best_time = -1;
    int K_hat = -1;
    int km_iter_best = 1;
    int km_iter_total = 1;
    int new_best_km = 0;

    // Dual Certificate
    int dc = 0;
    int dc_grad_iter = -1;
    double dc_time = -1;
    double* Y_a_best = (double *) R_alloc(d,sizeof(double));

    // Projected Gradient Descent A - declarations
    clock_t start_time = clock();
    clock_t cur_time;
    int grad_iter_best = -1;
    double grad_iter_best_time = 0;
    int grad_iter_total = 1;
    int outer_iterations = 0;
    int next_iter = 1;
    int num_momentum_restarts = 0;
    double lambda_min_tp1;
    double lambda_min_best;
    double* X_t = NULL;
    double* X_tp1 = NULL;
    double* Z_t = NULL;
    double* Z_tp1 = NULL; // use result as workspace
    double* Z_best = results->Z_best; // store actual result
    double lambda_t = 0; // auxiliary sequence for nesterov acc
    double lambda_tp1 = 1;
    double gamma_t;
    double obj_best;
    double obj_tp1;
    double obj_t;
    int stop_loop = 0;
    int early_stop = 0;
    int number_restarts = opts->number_restarts;
    int* restarts = opts->restarts;
    int current_restart = -1;
    int next_restart = -1;
    int last_restart = 0;
    double* GX_t = NULL;
    double* GS_t = NULL;
    double alpha = opts->alpha;

    if(number_restarts > 0){
        next_restart = *restarts;
        restarts++;
        current_restart = 1;
    }

    // Options
    int km_rep = opts->kmeans_iter;
    int verbosity = opts->verbosity;
    int finish_pgd = opts->finish_pgd;
    int primal_only = opts -> primal_only;
    int max_iter = opts -> max_iter;
    int dual_frequency = opts->dual_frequency;
    if(primal_only != 0) dual_frequency = max_iter + 1;
    int early_stop_mode = opts -> early_stop_mode;
    int early_stop_lag = opts -> early_stop_lag;
    double early_stop_eps = opts -> early_stop_eps;
    double* last_es_obj;
    if(early_stop_mode > 0){
        last_es_obj = (double *) R_alloc(early_stop_lag,sizeof(double));
        for(int i = 0; i < early_stop_lag; i++) {
            last_es_obj[i] = 0;
        }
    }

    ////////////////////////////////////////////////////////////
    //// STEP 2 - Initial K-means Solution, Certificate
    ////////////////////////////////////////////////////////////
    if(primal_only == 0) {
        tmp_hc_sol.clusters = km_clusters_best;
        hclust_FORCE(D_kmeans,d,&tmp_hc_sol,&work);
        K_hat = tmp_hc_sol.K;
        km_val_best = clust_to_opt_val_adapt(&prob,&tmp_hc_sol,&work);
        cur_time = clock();
        km_best_time = time_difference_ms(start_time,cur_time);
        kmeans_dual_solution_nok_impl(km_clusters_best,&prob,DUAL_EPS1_DEFAULT,Y_a_best, &dc, &work);
        if(dc == 1){
            cur_time = clock();
            dc_time = time_difference_ms(start_time,cur_time);
            dc_grad_iter = grad_iter_total;
        }
    }
    ////////////////////////////////////////////////////////////
    //// STEP 3 - Outer PGD Loop
    ////////////////////////////////////////////////////////////
    stop_loop = ((dc && !finish_pgd) || (grad_iter_total >= max_iter) );
    Z_tp1 = (double *) mem_pool_remove(&free_d2);
    X_tp1 = (double *) mem_pool_remove(&free_d2);
    memcpy(X_tp1,X0,d2*sizeof(double));
    memcpy(Z_tp1,X_tp1,d2*sizeof(double));
    smoothed_objective(&prob,Z_tp1,&lambda_min_tp1,&obj_tp1,&work);
    lambda_min_best = lambda_min_tp1;
    obj_best = obj_tp1;

    if(verbosity > -1){
        Rprintf("\tSolving K-Means SDP with FORCE\r\n");
        Rprintf("\t\tOptions -- Verbosity: %d\r\n",verbosity);
        Rprintf("\t\tOptions -- Finish PGD: %d\r\n",finish_pgd);
        Rprintf("\t\tOptions -- Max. Iter.: %d\r\n",max_iter);
        Rprintf("\t\tOptions -- Dual Frequency: %d\r\n",dual_frequency);
        Rprintf("\t\tOptions -- K-means Reps.: %d\r\n",km_rep);
        Rprintf("\t\tOptions -- Eps. Objective: %.3f\r\n",eps_obj);
        Rprintf("\t\tOptions -- Alpha: %.5f\r\n",alpha);
    }

    if(verbosity > 0){
        Rprintf("\tINITIAL VALUES - lambda_min_E(X_1) = %4.4f\r\n",lambda_min_tp1);
        Rprintf("\tINITIAL VALUES - f_mu,E(X_1) = %4.4f\r\n",obj_tp1);
    }

    while(!stop_loop){
        // STEP 3A -- Inner Loop
        next_iter = 1;
        outer_iterations++;
        if(verbosity > 0){
            Rprintf("\tOUTER ITERATION %d -- START\r\n",outer_iterations);
        }
        while(((grad_iter_total % dual_frequency > 0) || next_iter) && early_stop == 0) {
            next_iter = 0;
            if(verbosity > 1){
                Rprintf("\t\tINNER ITERATION %d -- START\r\n",grad_iter_total);
            }

            // STEP 3AA -- Check if time to restart X_t
            if(grad_iter_total == next_restart){
                if(verbosity > 1){
                    Rprintf("\t\tINNER ITERATION %d -- RESTART X0\r\n",grad_iter_total);
                }
                // No longer need X_t, X_tp1, Z_t, Z_tp1
                // X_tp1 <- E + ( 1/(1 - lambda_min_best) ) * (Z_best - E)
                project_E(&prob,Z_best,lambda_min_best,X_tp1);
                //Z_tp1 <- X_tp1
                memcpy(Z_tp1,X_tp1,d2*sizeof(double));
                // lambda_t <- 0
                lambda_t = 0;
                // lambda_tp1 <- 1
                lambda_tp1 = 1;
                // obj_tp1 <- s_res$objective_value
                // lambda_min_tp1 <- s_res$lambda_min
                smoothed_objective(&prob,Z_tp1,&lambda_min_tp1,&obj_tp1,&work);
                //set best result
                memcpy(Z_best,Z_tp1,d2*sizeof(double));
                obj_best = obj_tp1;
                lambda_min_best = lambda_min_tp1;
                grad_iter_best = grad_iter_total;
                last_restart = grad_iter_total;

                if(current_restart < number_restarts){
                    current_restart++;
                    next_restart = *restarts;
                    restarts++;
                }
            }

            // STEP 3AB - Update History
            X_t = X_tp1;
            Z_t = Z_tp1;
            X_tp1 = 0; //Can remove later
            Z_tp1 = 0; //Can remove later
            obj_t = obj_tp1;

            // STEP 3AC -- Update Auxiliary Sequences
            lambda_t = lambda_tp1;
            lambda_tp1 = (1+ sqrt(1 + 4*lambda_t*lambda_t))/2;
            gamma_t = (1-lambda_t) / lambda_tp1;

            // STEP 3AD -- Find Gradient
            GX_t = (double *) mem_pool_remove(&free_d2);
            GS_t = (double *) mem_pool_remove(&free_d2);
            smoothed_gradient(&prob, X_t, GX_t, GS_t, &work);

            // STEP 3AE -- Update Primary Sequences
            C_perp_update_nok(&prob,alpha,X_t,GX_t,GS_t,&work);
            Z_tp1 = X_t;
            X_t = 0;
            smoothed_objective(&prob,Z_tp1,&lambda_min_tp1,&obj_tp1,&work);

            // Memory @ GS_t can be freed
            // Memory @ X_t can be freed
            mem_pool_insert(&free_d2,GS_t);
            mem_pool_insert(&free_d2,GX_t);
            GX_t = 0;
            GS_t = 0;

            // Update X_tp1
            //  X_tp1 -> (1 - gamma_t)*Z_tp1 + gamma_t*Z_t
            F77_NAME(dscal)(&d2,&gamma_t,Z_t,&INC1);
            dtmp1 = 1 - gamma_t;
            F77_NAME(daxpy)(&d2,&dtmp1,Z_tp1,&INC1,Z_t,&INC1);
            X_tp1 = Z_t;
            Z_t = 0;

            if(verbosity > 1){
                Rprintf("\t\tINNER ITERATION %d -- Smoothed Objective: %4.4f\r\n",grad_iter_total,obj_tp1);
                Rprintf("\t\tINNER ITERATION %d -- Lambda_min: %4.4f\r\n",grad_iter_total,lambda_min_tp1);
                Rprintf("\t\tINNER ITERATION %d -- lambda_t: %4.4f\r\n",grad_iter_total,lambda_t);
                Rprintf("\t\tINNER ITERATION %d -- lambda_tp1: %4.4f\r\n",grad_iter_total,lambda_tp1);
                Rprintf("\t\tINNER ITERATION %d -- gamma_t: %4.4f\r\n",grad_iter_total,gamma_t);
            }
            if(verbosity == 5){
                // dtmp1 = clust_to_opt_val(&prob,km_clusters_best,&work);
                dtmp1 = 0.0;
                Rprintf("\t\tINNER ITERATION %d -- KM SDP Objective: %4.4f\r\n",grad_iter_total,dtmp1);
            }

            // STEP 3AF -- Check for Adaptive Restart of Aux Sequences
            if(obj_tp1 < obj_t) {
                lambda_t = 0.0;
                lambda_tp1 = 1.0;
                num_momentum_restarts++;
            }

            // STEP 3AG -- Update Best Solution
            if(obj_tp1 > obj_best){
                if(verbosity > 1){
                    Rprintf("\t\tINNER ITERATION %d -- Found New Best\r\n",grad_iter_total);
                }
                memcpy(Z_best,Z_tp1,d2*sizeof(double));
                obj_best = obj_tp1;
                lambda_min_best = lambda_min_tp1;
                grad_iter_best = grad_iter_total;
                cur_time = clock();
                grad_iter_best_time = time_difference_ms(start_time,cur_time);
            }

            if(verbosity > 0){
                Rprintf("\t\tINNER ITERATION %d -- COMPLETE\r\n",grad_iter_total);
            }
            grad_iter_total++;

            // STEP 3AH -- Check early stop relative error criterion
            if(early_stop_mode > 0) {
                last_es_obj[ grad_iter_total % early_stop_lag] = obj_tp1;
            }
            if(early_stop_mode == 1 && grad_iter_total-last_restart > early_stop_lag) {
                //absolute error
                dtmp1 = min_array(early_stop_lag,last_es_obj);
                dtmp2 = max_array(early_stop_lag,last_es_obj);
                dtmp1 = dtmp2 - dtmp1;
                early_stop = dtmp1 > early_stop_eps ? 0 : 1;
            }
            if(early_stop_mode == 2 && grad_iter_total-last_restart > early_stop_lag) {
                //relative error
                dtmp1 = min_array(early_stop_lag,last_es_obj);
                dtmp2 = max_array(early_stop_lag,last_es_obj);
                dtmp1 = (dtmp2 - dtmp1) / (dtmp1 + 0.000001);
                early_stop = dtmp1 > early_stop_eps ? 0 : 1;
            }

        }

        // STEP 3B -- Dual Certificate Search -- HC is deterministic
        new_best_km = 0;
        project_E(&prob,Z_best,lambda_min_best,results->B_Z_best);
        tmp_hc_sol.clusters = km_clusters_new;
        hclust_FORCE(results->B_Z_best,d,&tmp_hc_sol,&work);
        km_val_new = clust_to_opt_val_adapt(&prob,&tmp_hc_sol,&work);
        km_iter_total++;
        if(km_val_new < km_val_best) {
            km_val_best = km_val_new;
            km_clusters_tmp = km_clusters_best;
            km_clusters_best = km_clusters_new;
            km_clusters_new = km_clusters_tmp;
            km_iter_best = km_iter_total;
            K_hat = tmp_hc_sol.K;
            new_best_km = 1;
            cur_time = clock();
            km_best_time = time_difference_ms(start_time,cur_time);
        }
        if(new_best_km && dc == 0){
            kmeans_dual_solution_nok_impl(km_clusters_best,&prob,DUAL_EPS1_DEFAULT,Y_a_best, &dc, &work);
            if(dc == 1){
                cur_time = clock();
                dc_time = time_difference_ms(start_time,cur_time);
                dc_grad_iter = grad_iter_total;
            }
        }

        if(verbosity > 0){
            Rprintf("\tOUTER ITERATION %d -- Dual Feasible: %d\r\n",outer_iterations,dc);
            Rprintf("\tOUTER ITERATION %d -- COMPLETE\r\n",outer_iterations);
        }

        // STEP 3C - Update Stopping Criterion
        stop_loop = ((dc && !finish_pgd) || (grad_iter_total >= max_iter) || early_stop == 1);
    }

    ////////////////////////////////////////////////////////////
    //// STEP 4 - Set Return Values
    ////////////////////////////////////////////////////////////
    // Copy Z_T for output, and find projection
    memcpy(results->Z_T,Z_tp1,d2*sizeof(double));
    project_E(&prob,Z_tp1,lambda_min_tp1,results->B_Z_T);

    // Project Z_best
    project_E(&prob,Z_best,lambda_min_best,results->B_Z_best);

    //Find optimal values
    results->B_Z_T_opt_val = F77_NAME(ddot)(&d2,D,&INC1,results->B_Z_T,&INC1);
    results->B_Z_best_opt_val = F77_NAME(ddot)(&d2,D,&INC1,results->B_Z_best,&INC1);

    // Copy out kmeans best
    if(primal_only == 0) memcpy(results->kmeans_best,km_clusters_best,d*sizeof(int));

    // Set scalar return values
    results->Z_T_lmin = lambda_min_tp1;
    results->Z_best_lmin = lambda_min_best;
    results->kmeans_best_time = km_best_time;
    results->kmeans_iter_best = km_iter_best;
    results->kmeans_iter_total = km_iter_total;
    results->dc = dc;
    results->dc_time = dc_time;
    results->dc_grad_iter = dc_grad_iter;
    results->grad_iter_best = grad_iter_best;
    results->grad_iter_best_time = grad_iter_best_time;
    results->kmeans_opt_val = km_val_best;
    cur_time = clock();
    results->total_time = time_difference_ms(start_time,cur_time);

    if(verbosity > -1){
        Rprintf("\tFORCE Algorithm Complete\r\n");
        Rprintf("\t\tDual Certificate:%d\r\n",dc);
        Rprintf("\t\tInner Iterations Total:%d\r\n",grad_iter_total);
        Rprintf("\t\tOuter Iterations Total:%d\r\n",outer_iterations);
        Rprintf("\t\tTotal Running Time:%.3fs\r\n",results->total_time);
        Rprintf("\t\t<D,B_Z_T> = %.4f\r\n",results->B_Z_T_opt_val);
        Rprintf("\t\t<D,B_Z_best> = %.4f\r\n",results->B_Z_best_opt_val);
        if(primal_only == 0) Rprintf("\t\t<D,B_km> = %.4f\r\n",km_val_best);
    }
    if(verbosity == 5){
        Rprintf("\t\tTotal Momentum Restarts:%d\r\n",num_momentum_restarts);
    }
}


// Access point for R code --- needed to pass in options because
// cannot pass struct
void FORCE_adapt_R(double* D, double* D_kmeans, double* E, double* ESI, double* X0, int* d,
    int* in_verbosity, int* in_kmeans_iter, int* in_dual_frequency, int* in_max_iter,
    int* in_finish_pgd, int* in_primal_only, int* in_number_restarts, int* in_restarts, double* in_alpha, double* in_eps_obj,
    int* in_early_stop_mode, int* in_early_stop_lag, double* in_early_stop_eps,
    double* out_Z_T, double* out_B_Z_T, double* out_Z_T_lmin, double* out_Z_best, double* out_B_Z_best, double* out_Z_best_lmin,
    double* out_B_Z_T_opt_val, double* out_B_Z_best_opt_val, double* out_kmeans_opt_val,  int* out_kmeans_best, double* out_kmeans_best_time,
    int* out_kmeans_iter_best, int* out_kmeans_iter_total, int* out_dc, double* out_dc_time,
    int* out_dc_grad_iter, int* out_grad_iter_best, double* out_grad_iter_best_time, double* out_total_time)
{
    pgd_opts opts_in;
    pgd_results results;

    // populate options structure, and return pointers
    opts_in.verbosity = *in_verbosity;
    opts_in.kmeans_iter = *in_kmeans_iter;
    opts_in.dual_frequency = *in_dual_frequency;
    opts_in.max_iter = *in_max_iter;
    opts_in.finish_pgd = *in_finish_pgd;
    opts_in.number_restarts = *in_number_restarts;
    opts_in.restarts = in_restarts;
    opts_in.alpha = *in_alpha;
    opts_in.eps_obj = *in_eps_obj;
    opts_in.primal_only = *in_primal_only;
    opts_in.early_stop_mode = *in_early_stop_mode;
    opts_in.early_stop_lag = *in_early_stop_lag;
    opts_in.early_stop_eps = *in_early_stop_eps;

    results.Z_T = out_Z_T;
    results.B_Z_T = out_B_Z_T;
    results.Z_best = out_Z_best;
    results.B_Z_best = out_B_Z_best;
    results.kmeans_best = out_kmeans_best;

    // call primal_dual_adar
    FORCE_adapt(D,D_kmeans,E,ESI,X0,*d, &opts_in, &results);

    // assign out_* values from results
    // out_Z_T passed by reference
    // out_Z_best passed by reference
    // out_kmeans_best passed by reference
    *out_Z_T_lmin = results.Z_T_lmin;
    *out_Z_best_lmin = results.Z_best_lmin;
    *out_kmeans_best_time = results.kmeans_best_time;
    *out_kmeans_iter_best = results.kmeans_iter_best;
    *out_kmeans_iter_total = results.kmeans_iter_total;
    *out_dc = results.dc;
    *out_dc_time = results.dc_time;
    *out_dc_grad_iter = results.dc_grad_iter;
    *out_grad_iter_best = results.grad_iter_best;
    *out_grad_iter_best_time = results.grad_iter_best_time;
    *out_kmeans_opt_val = results.kmeans_opt_val;
    *out_B_Z_T_opt_val = results.B_Z_T_opt_val;
    *out_B_Z_best_opt_val = results.B_Z_best_opt_val;
    *out_total_time = results.total_time;
}
