#ifndef UTIL_H
#define UTIL_H

#include "time.h"
#include "FORCE.h"

// CONSTANTS
extern const double DUAL_EPS1_DEFAULT, DUAL_EPS2_DEFAULT, DUAL_Y_T_MIN_DEFAULT,GNORM_FACTOR,ALPHA,BETA;
extern const int INC1,NRHS1;
extern const char JOBZV,JOBZN,UPLO,SIDE_L,SIDE_R,TRANS_N,TRANS_T;

// NEW TYPES
typedef struct mem_pool {
    void** base;
    int length;
    int start_idx;
    int end_idx;
} mem_pool;

typedef struct workspace {
    double* dwork;
    int* iwork;
    int ldwork;
    int liwork;
    int dsyevd_liwork;
    int dsyevd_liwork_N;
    int dsyevd_ldwork;
    int dsyevd_ldwork_N;

} workspace;


typedef struct problem_instance {
    double* D;
    double* E;
    double* ESI;
    double* D_rsums;
    double TrD;
    double DTD;
    double mu_n;
    int d;
    int K;
} problem_instance;


// External Functions

// GENERAL UTILITIES, MEMORY MANAGEMENT
double min_array(int d,double* V);
double max_array(int n, double* V);
void* mem_pool_remove(mem_pool* pool);
void mem_pool_insert(mem_pool* pool, void* mem_ptr);
void allocate_workspace_FORCE(int d, int K, workspace* work);
void allocate_workspace_FORCE_adapt(int d, workspace* work);
void initialize_problem_instance(double* D, double* E, double* ESI, double mu,
                                int d,int K, problem_instance* prob);
double time_difference_ms(clock_t start, clock_t end);
void initialize_identity_matrix(double* restrict I, int d);

// KMEANS CLUSTERING
void kmeans_pp_impl(double const* D, int K, int n, int m, int* cluster_assignment_r,
                        double* centers_r, int* num_iters_R, double* time_R, workspace* work);

// HC CLUSTERING


#endif
