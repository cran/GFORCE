#include "float.h"
#include "math.h"
#include "time.h"
#include "R.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "util.h"
#include "util_mops.h"

const double DUAL_EPS1_DEFAULT = 0.01;
const double DUAL_EPS2_DEFAULT = 0.000001;
const double DUAL_Y_T_MIN_DEFAULT = 0.01;
const char JOBZV = 'V';
const char JOBZN = 'N';
const char UPLO = 'L';
const char SIDE_L = 'L';
const char SIDE_R = 'R';
const char TRANS_N = 'N';
const char TRANS_T = 'T';
const double ALPHA = 1.0;
const double BETA = 0.0;
const int INC1 = 1;
const int NRHS1 = 1;
const double GNORM_FACTOR = 2.82842712475; //\sqrt{8}

// GENERAL UTILITIES, MEMORY MANAGEMENT
double min_array(int n, double* V){
    double min_val = DBL_MAX;
    for(int i=0; i < n; i++) min_val = (V[i] < min_val) ? V[i] : min_val;
    return min_val;
}

double max_array(int n, double* V){
    double max_val = -DBL_MAX;
    for(int i=0; i < n; i++) max_val = (V[i] > max_val) ? V[i] : max_val;
    return max_val;
}


void mem_pool_insert(mem_pool* pool, void* mem_ptr){
    int tmp_idx = (pool->end_idx);
    if(tmp_idx == (pool->length) -1) {
       tmp_idx = -1;
    }
    *((pool->base) + ++tmp_idx) = mem_ptr;
    pool->end_idx = tmp_idx;
}

void* mem_pool_remove(mem_pool* pool){
    int tmp_idx = pool->start_idx;
    void* mem_ptr = *((pool->base) + tmp_idx++);
    if(tmp_idx == pool->length){
        tmp_idx = 0;
    }
    pool->start_idx = tmp_idx;
    return mem_ptr;
}





// Allocates properly sized workspace for primal_dual_adar
void allocate_workspace_FORCE(int d, int K, workspace* work){
    // local vars
    int liwork = -1;
    int liwork_N = -1;
    int ldwork = -1;
    int ldwork_N = -1;
    int tmp1;

    // USE CASE 0 -- Find Space Needed for both types of DSYEVD Calls
    double* X = NULL;
    double dwork0[2];
    int iwork0[2];
    double dvec[2];
    int lapack_info = 0;
    F77_CALL(dsyevd)(&JOBZN,&UPLO,&d,X,&d,dvec,dwork0,&ldwork_N,iwork0,&liwork_N,&lapack_info);
    liwork_N = *iwork0;
    ldwork_N = (int) *dwork0;
    work -> dsyevd_liwork_N = liwork_N;
    work -> dsyevd_ldwork_N = ldwork_N;
    F77_CALL(dsyevd)(&JOBZV,&UPLO,&d,X,&d,dvec,dwork0,&ldwork,iwork0,&liwork,&lapack_info);
    liwork = *iwork0;
    ldwork = (int) *dwork0;
    work -> dsyevd_liwork = liwork;
    work -> dsyevd_ldwork = ldwork;


    // USE CASE 1 -- SMOOTHED_GRADIENT COMPUTATION
    ldwork = ldwork + d; //add space for w to hold eigvals

    // Workspace size for smoothed gradient (ignoring) syevd is 2d^2+d
    tmp1 = 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;


    // USE CASE 2 -- C_perp_update and clust_to_opt_val
    // Get workspace size needed for C_perp_update and clust_to_opt_val
    // These require ldwork >= d^2+5d+2 and liwork >= d+3K+3
    tmp1 = d*d + 5*d + 2;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    tmp1 = d + 3*K + 3;
    liwork = liwork > tmp1 ? liwork : tmp1;
    

    // USE CASE 3 -- kmeans_pp
    // STRICTLY DOMINATED IN REQUIREMENTS BY CASE 2


    // USE CASE 4 -- smoothed_objective
    // STRICTLY DOMINATED BY REQUIREMENTS IN USE CASE 2 (should be)
    tmp1 = ldwork_N + 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    liwork = liwork > liwork_N ? liwork : liwork_N;


    // Initialize workspace
    work -> dwork = (double *) R_alloc(ldwork,sizeof(double));
    work -> ldwork = ldwork;
    work -> iwork = (int *) R_alloc(liwork,sizeof(int));
    work -> liwork = liwork;
}

// Allocates properly sized workspace for primal_dual_adar
void allocate_workspace_FORCE_adapt(int d, workspace* work){
    // local vars
    int liwork = -1;
    int liwork_N = -1;
    int ldwork = -1;
    int ldwork_N = -1;
    int tmp1;

    // USE CASE 0 -- Find Space Needed for both types of DSYEVD Calls
    double* X = NULL;
    double dwork0[2];
    int iwork0[2];
    double dvec[2];
    int lapack_info = 0;
    F77_CALL(dsyevd)(&JOBZN,&UPLO,&d,X,&d,dvec,dwork0,&ldwork_N,iwork0,&liwork_N,&lapack_info);
    liwork_N = *iwork0;
    ldwork_N = (int) *dwork0;
    work -> dsyevd_liwork_N = liwork_N;
    work -> dsyevd_ldwork_N = ldwork_N;
    F77_CALL(dsyevd)(&JOBZV,&UPLO,&d,X,&d,dvec,dwork0,&ldwork,iwork0,&liwork,&lapack_info);
    liwork = *iwork0;
    ldwork = (int) *dwork0;
    work -> dsyevd_liwork = liwork;
    work -> dsyevd_ldwork = ldwork;


    // USE CASE 1 -- SMOOTHED_GRADIENT COMPUTATION
    ldwork = ldwork + d; //add space for w to hold eigvals

    // Workspace size for smoothed gradient (ignoring) syevd is 2d^2+d
    tmp1 = 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;


    // USE CASE 2 -- C_perp_update and clust_to_opt_val
    // Get workspace size needed for C_perp_update and clust_to_opt_val
    // These require ldwork >= d^2+5d+2 and liwork >= d+3K+3 (K<= d)
    tmp1 = d*d + 5*d + 2;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    tmp1 = d + 3*d + 3;
    liwork = liwork > tmp1 ? liwork : tmp1;
    

    // USE CASE 3 -- Hclust
    // ldwork >= d^2 + 2d, liwork >= 7d+3
    tmp1 = 3*d*d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    tmp1 = 7*d + 3;
    liwork = liwork > tmp1 ? liwork : tmp1;

    // USE CASE 4 -- smoothed_objective
    // STRICTLY DOMINATED BY REQUIREMENTS IN USE CASE 2 (should be)
    tmp1 = ldwork_N + 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    liwork = liwork > liwork_N ? liwork : liwork_N;


    // Initialize workspace
    work -> dwork = (double *) R_alloc(ldwork,sizeof(double));
    work -> ldwork = ldwork;
    work -> iwork = (int *) R_alloc(liwork,sizeof(int));
    work -> liwork = liwork;
}


// Takes input and precomputes values for problem instance
// computes the fields D_rsums, TrD, DTD
// Assumes that these are un-initialized
void initialize_problem_instance(double* D, double* E, double* ESI, double mu, int d,
                                    int K, problem_instance* prob){
    int d2 = d*d;
    double* D_rsums = (double *) R_alloc(d,sizeof(double));
    dcsum(D,d,D_rsums);

    prob -> d = d;
    prob -> D = D;
    prob -> E = E;
    prob -> ESI = ESI;
    prob -> D_rsums = D_rsums;
    prob -> mu_n = -1*mu;
    prob -> K = K;
    prob -> TrD = dtrace(D,d);
    prob -> DTD = F77_NAME(ddot)(&d2,D,&INC1,D,&INC1);
}

void initialize_identity_matrix(double* restrict I,int d){
    int d2 = d*d;
    int dp1 = d+1;
    double* ptmp1 = I;
    for (int i=0; i < d2; i++){
        *ptmp1 = 0;
        ptmp1++;
    }
    ptmp1 = I;
    for(int i=0; i < d; i++){
        *ptmp1 = 1;
        ptmp1 = ptmp1 + dp1;
    }
}


double time_difference_ms(struct timespec* start,struct timespec* end){
    double diff;
    diff = (end->tv_sec - start->tv_sec);
    diff += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return diff;
}


double euclidean_distance(double const* restrict p1, double const* restrict p2, int m, double* restrict euclidean_distance_tmp) {
    double acc = 0.0;

    #pragma omp simd
    for(int i=0; i < m; i++){
        //it seems to vectorize better this way...
        euclidean_distance_tmp[i] = p1[i] - p2[i];
        acc = acc + euclidean_distance_tmp[i] * euclidean_distance_tmp[i] ;
    } 

    return acc;
}
