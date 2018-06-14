#ifndef UTIL_FORCE_H
#define UTIL_FORCE_H
#include "util.h"
#include "FORCE.h"

// External Functions

// FORCE PRIMAL METHOD HELPERS
double clust_to_opt_val(problem_instance* prob, int* ga_hat, workspace* work);
double clust_to_opt_val_adapt(problem_instance* prob,hclust_t* tmp_hc_sol,workspace* work);
void smoothed_gradient(problem_instance* prob, double* X, double* GX_t, double* GS_t,
                        workspace* work);
void smoothed_gradient_nok(problem_instance* prob, double* X, double* GX_t, double* GS_t,
                        workspace* work);
void smoothed_objective(problem_instance* prob, double* X, double* lambda_min,
                        double* obj_val, workspace* work);
void smoothed_objective_nok(problem_instance* prob, double* X, double* lambda_min,
                        double* obj_val, workspace* work);
void C_perp_update(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace* work);
void C_perp_update_nok(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace*work);
void project_E(problem_instance* prob, double* Z, double lmin, double* Z_proj);


// DUAL CERTIFICATES
void kmeans_dual_solution_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1, double eps2,
                                double Y_T_min, double* restrict Y_a_r, double* restrict Y_T_r, int* restrict feasible_r,
                                workspace* restrict work);
void kmeans_dual_solution_nok_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1,
                                   double* restrict Y_a_r, int* restrict feasible_r, workspace* restrict work);

void hclust_FORCE(double* dists,int d,hclust_t* hclust_sol,workspace* work);

#endif
