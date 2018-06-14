#ifndef TEST_HOOKS_H
#define TEST_HOOKS_H
#include "util.h"
// #include "util_FORCE.h"

//////////////////////////////////////////////////////////////////
////// UNIT TEST R ENTRY POINTS FOR LOW LEVEL & UTILITY FUNCTIONS
//////////////////////////////////////////////////////////////////

// primal_dual_adar_helpers.c
void test_smoothed_gradient_S_base(double* X, double* E,double* GS_t, int* d, double* S_min_r);
void test_smoothed_gradient_X_base(double* X, double* ESI, double* GX_t, int* d, int* K, double* X_min_r);
void test_smoothed_gradient(double* X, double* E,double* ESI,int* d0, int* K0, double* mu0,double* GX_t,double* GS_t);
void test_project_C_perpendicular(double* D, int* d0, int* K0, double* GX_t, double* GS_t);
void test_project_C_perpendicular_nok(double* D, int* d0, double* GX_t, double* GS_t);
void test_project_E(double* E, double* Z, int* d, double* lmin, double* Z_proj);
void test_smoothed_objective(double* X, double* E,double* ESI,int* d0, int* K0,
                            double* mu0, double* lambda_min, double* obj_val);
void test_clust_to_opt_val(double* D, int* d, int* K, int* clusters, double* opt_val_r);


// convex_kmeans_util.c
void test_daps(double* A,double* c, int* d);
void test_dsmtd(double* A,double* B, int* d, char** side);
void test_dvexp(double* A, int* d);
void test_dsumv(double* A, int* d, double* sum_r);
void test_dtrace(double* A, int* d, double* trace_r);
void test_dcsum(double* A, int* d, double* A_csums);
void test_dxpyez(int* d, double* X, double* Y, double* Z);


/////////////////////////////////////////////////////////////////////////
////// UNIT TEST C ENTRY POINTS FOR LOW LEVEL FUNCTIONS W/OUT ENTRYPOINT
/////////////////////////////////////////////////////////////////////////

void smoothed_gradient_S_base(double* X, double* E,double* GS_t, int d, double* S_min_r);
void smoothed_gradient_X_base(double* X, double* ESI, double* GX_t, double* d2_tmp,
                                workspace* work, int d, int X_eigs_idx);
void project_C_perpendicular(problem_instance* prob, double* GX_t, double* GS_t,
                                workspace* work);
void project_C_perpendicular_nok(problem_instance* prob, double* GX_t, double* GS_t,
                                  workspace* work);
#endif
