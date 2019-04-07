#include "math.h"
#include "stdio.h"
#include "float.h"
#include "R.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"
#include "util_FORCE.h"
#include "util_mops.h"
#include "test_hooks.h"

void smoothed_objective_X_base(double* X, double* ESI, double* d2_tmp, double* d2_tmp2,
                                workspace* work, int d, int X_eigs_idx);
void smoothed_gradient_X_base(double* X, double* ESI, double* GX_t, double* d2_tmp,
                                workspace* work, int d, int X_eigs_idx);
void smoothed_gradient_S_base(double* X, double* E,double* GS_t, int d, double* S_min_r);


// NOTE: All methods assume column major layout

// CANNOT OVERWRITE *X or *E or *ESI
// REQUIRES work -> dwork be of length at least min(dseyvd_ldwork + d,
//   2d^2 + d) (total should be 1 + 7d +2d^2)
// REQUIRES work -> iwork be of length at least dseyvd_ldwork (should be 3+5d)
void smoothed_gradient(problem_instance* prob, double* X, double* GX_t, double* GS_t, 
                        workspace* work){

    // Grab Problem Parameters
    double* E = prob -> E;
    double* ESI = prob -> ESI;
    int d = prob -> d;
    double mu_n = prob -> mu_n;

    // Local Vars
    double* d2_tmp;
    double* d2_tmp2;
    double S_min;
    double X_min;
    double lambda_min_n;
    double scale_factor;
    double* X_eigs;
    int X_eigs_idx;
    int tmp;
    int d2 = d*d;
    
    X_eigs_idx = work -> dsyevd_ldwork;;
    tmp = 2*d2;
    X_eigs_idx = X_eigs_idx > tmp ? X_eigs_idx : tmp;
    X_eigs = (work -> dwork) + X_eigs_idx;
    d2_tmp = work -> dwork;
    d2_tmp2 = d2_tmp + d2;

    smoothed_gradient_X_base(X,ESI,GX_t,GS_t,work,d,X_eigs_idx);
    smoothed_gradient_S_base(X,E,GS_t,d,&S_min);

    X_min = min_array(d,X_eigs);
    lambda_min_n = -1*(S_min < X_min ? S_min : X_min);

    // subtract lambda_min from all eigenvalues
    daps(X_eigs,1,lambda_min_n,d);
    daps(GS_t,1,lambda_min_n,d2);

    // Divide both eigenvalues by mu
    F77_NAME(drscl)(&d,&mu_n,X_eigs,&INC1);
    F77_NAME(drscl)(&d2,&mu_n,GS_t,&INC1);

    // Exponentiate eigenvalues
    dvexp(X_eigs,d);
    dvexp(GS_t,d2);

    // Compute Scaling Factor
    scale_factor = dsumv(X_eigs,d);
    scale_factor = scale_factor + dsumv(GS_t,d2);

    // Rescale X eigenvalues
    memcpy(d2_tmp,GX_t,d2*sizeof(double));
    memcpy(d2_tmp2,GX_t,d2*sizeof(double));
    dsmtd(d2_tmp,X_eigs,d,SIDE_R);//multiply V by diag(X_eigs), stored in dtmp2
    //neither matrix is symmetric, so need to use dgemm
    F77_NAME(dgemm)(&TRANS_N,&TRANS_T,&d,&d,&d,&ALPHA,d2_tmp,&d,d2_tmp2,&d,&BETA,GX_t,&d);

    // Rescale all eigenvalues by scaling factor
    F77_NAME(drscl)(&d2,&scale_factor,GX_t,&INC1);
    F77_NAME(drscl)(&d2,&scale_factor,GS_t,&INC1);
}

// X is untouched
// REQUIRES dwork 2d^2 + 3d + 1 (assuming dsyevd only needs 2d+1)
// REQUIRES iwork 1.
void smoothed_objective(problem_instance* prob, double* X, double* lambda_min_tp1,
                        double* obj_tp1, workspace* work){
    // Grab Problem Parameters
    double* E = prob -> E;
    double* ESI = prob -> ESI;
    int d = prob -> d;
    double mu_n = prob -> mu_n;

    // Local Vars
    double* GS_t;
    double* d2_tmp;
    double S_min;
    double X_min;
    double lambda_min_n;
    double lambda_min;
    double obj_value;
    double* X_eigs;
    int X_eigs_idx;
    int d2 = d*d;
    
    X_eigs_idx = work -> dsyevd_ldwork_N;
    X_eigs = (work -> dwork) + X_eigs_idx;
    GS_t = X_eigs + d;
    d2_tmp = GS_t + d2;

    smoothed_objective_X_base(X,ESI,GS_t,d2_tmp,work,d,X_eigs_idx);
    smoothed_gradient_S_base(X,E,GS_t,d,&S_min);

    // Get lambda_min and negate
    X_min = min_array(d,X_eigs);
    lambda_min = S_min < X_min ? S_min : X_min;
    lambda_min_n = -1*lambda_min;
    *lambda_min_tp1 = lambda_min;

    // subtract lambda_min from all eigenvalues
    daps(X_eigs,1,lambda_min_n,d);
    daps(GS_t,1,lambda_min_n,d2);

    // Divide both eigenvalues by mu
    F77_NAME(drscl)(&d,&mu_n,X_eigs,&INC1);
    F77_NAME(drscl)(&d2,&mu_n,GS_t,&INC1);

    // Exponentiate eigenvalues
    dvexp(X_eigs,d);
    dvexp(GS_t,d2);

    // Compute Objective Value
    obj_value = dsumv(X_eigs,d);
    obj_value = obj_value + dsumv(GS_t,d2);
    obj_value = mu_n*log(obj_value) + lambda_min;

    *obj_tp1 = obj_value;
}


void smoothed_gradient_S_base(double* X, double* E,double* GS_t, int d, double* S_min_r){
    double S_min = DBL_MAX;
    double tmp;
    for(int i=0; i < d*d; i++){
        tmp = *X / *E;
        S_min = S_min > tmp ? tmp : S_min;
        *GS_t = tmp;
        X++;
        E++;
        GS_t++;
    }
    *S_min_r = S_min;
}

// REQUIRES work -> dwork be of length at least dseyvd_ldwork
//                  (should be 1 + 7d +2d^2) <- includes d entries for w
// REQUIRES work -> iwork be of length at least dseyvd_ldwork (should be 3+5d)
// Output is GX_t and index into work -> dwork where the eigenvalues can be found
void smoothed_gradient_X_base(double* X, double* ESI, double* GX_t, double* d2_tmp,
                                workspace* work, int d, int X_eigs_idx){
    double* w = (work->dwork) + X_eigs_idx;

    int lapack_info;
    F77_NAME(dsymm)(&SIDE_L,&UPLO,&d,&d,&ALPHA,X,&d,ESI,&d,&BETA,d2_tmp,&d);
    F77_NAME(dsymm)(&SIDE_L,&UPLO,&d,&d,&ALPHA,ESI,&d,d2_tmp,&d,&BETA,GX_t,&d); // computes ESI * X * ESI

    //eigenvectors will overwrite input matrix, so GX_t becomes V, eigenvectors
    F77_CALL(dsyevd)(&JOBZV,&UPLO,&d,GX_t,&d,w,work->dwork,&(work->dsyevd_ldwork),
                work->iwork,&(work->dsyevd_liwork),&lapack_info);
}

// REQUIRES work -> dwork be of length at least dseyvd_ldwork_N
//          (should be 1 + 2d +2d^2) <- includes d entries for w
// REQUIRES work -> iwork be of length at least dseyvd_ldwork_N (should be 1)
// output is w
void smoothed_objective_X_base(double* X, double* ESI, double* d2_tmp, double* d2_tmp2,
                                workspace* work, int d, int X_eigs_idx){
    double* w = (work->dwork) + X_eigs_idx;

    int lapack_info;
    F77_NAME(dsymm)(&SIDE_L,&UPLO,&d,&d,&ALPHA,X,&d,ESI,&d,&BETA,d2_tmp2,&d);
    F77_NAME(dsymm)(&SIDE_L,&UPLO,&d,&d,&ALPHA,ESI,&d,d2_tmp2,&d,&BETA,d2_tmp,&d); // computes ESI * X * ESI

    F77_CALL(dsyevd)(&JOBZN,&UPLO,&d,d2_tmp,&d,w,work->dwork,&(work->dsyevd_ldwork_N),
                work->iwork,&(work->dsyevd_liwork_N),&lapack_info);
}


// REQUIRES work -> dwork be of length K
// REQUIRES work -> iwork be of length d+3K+3
// supports clusters numbered 1..K or 0..K-1
double clust_to_opt_val(problem_instance* prob, int* ga_hat, workspace* work) {
    // Local Vars
    double opt_val;
    double* D = prob -> D;
    int d = prob -> d;
    int K = prob -> K;
    double* dwork = work -> dwork;
    int* iwork = work -> iwork;

    opt_val = dabgtp(D,ga_hat,d,K,iwork,dwork);

    return opt_val;
}

double clust_to_opt_val_adapt(problem_instance* prob,hclust_t* tmp_hc_sol,workspace* work) {
    // Local Vars
    double opt_val;
    double* D = prob -> D;
    int d = prob -> d;
    int K = tmp_hc_sol -> K;
    int* ga_hat = tmp_hc_sol -> clusters;
    double* dwork = work -> dwork;
    int* iwork = work -> iwork;

    opt_val = dabgtp(D,ga_hat,d,K,iwork,dwork);

    return opt_val;
}


// CANNOT OVERWRITE *X
// REQUIRES work -> dwork be of length at least min(dseyvd_ldwork + d,
//   2d^2 + d) (total should be 1 + 7d +2d^2)
// REQUIRES work -> iwork be of length at least dseyvd_ldwork (should be 3+5d)
void smoothed_gradient_nok(problem_instance* prob, double* X, double* GX_t, double* GS_t, 
                        workspace* work){

    // Grab Problem Parameters
    int d = prob -> d;
    double mu_n = prob -> mu_n;

    // Local Vars
    double* d2_tmp;
    double* d2_tmp2;
    double S_min;
    double X_min;
    double lambda_min_n;
    double scale_factor;
    double* X_eigs;
    int X_eigs_idx;
    int tmp;
    int d2 = d*d;
    int lapack_info;
    
    X_eigs_idx = work -> dsyevd_ldwork;
    tmp = 2*d2;
    X_eigs_idx = X_eigs_idx > tmp ? X_eigs_idx : tmp;
    X_eigs = (work -> dwork) + X_eigs_idx;
    d2_tmp = work -> dwork;
    d2_tmp2 = d2_tmp + d2;

    //eigenvectors will overwrite input matrix, so GX_t becomes V, eigenvectors
    memcpy(GX_t,X,d2*sizeof(double));
    F77_CALL(dsyevd)(&JOBZV,&UPLO,&d,GX_t,&d,X_eigs,work->dwork,&(work->ldwork),
                work->iwork,&(work->liwork),&lapack_info);
    // GS_t is just current value of X
    memcpy(GS_t,X,d2*sizeof(double));

    X_min = min_array(d,X_eigs);
    S_min = min_array(d2,GS_t);
    lambda_min_n = -1*(S_min < X_min ? S_min : X_min);

    // subtract lambda_min from all eigenvalues
    daps(X_eigs,1,lambda_min_n,d);
    daps(GS_t,1,lambda_min_n,d2);

    // Divide both eigenvalues by mu
    F77_NAME(drscl)(&d,&mu_n,X_eigs,&INC1);
    F77_NAME(drscl)(&d2,&mu_n,GS_t,&INC1);

    // Exponentiate eigenvalues
    dvexp(X_eigs,d);
    dvexp(GS_t,d2);

    // Compute Scaling Factor
    scale_factor = dsumv(X_eigs,d);
    scale_factor = scale_factor + dsumv(GS_t,d2);

    // Rescale X eigenvalues
    memcpy(d2_tmp,GX_t,d2*sizeof(double));
    memcpy(d2_tmp2,GX_t,d2*sizeof(double));
    dsmtd(d2_tmp,X_eigs,d,SIDE_R);//multiply V by diag(X_eigs), stored in dtmp2
    //neither matrix is symmetric, so need to use dgemm
    F77_NAME(dgemm)(&TRANS_N,&TRANS_T,&d,&d,&d,&ALPHA,d2_tmp,&d,d2_tmp2,&d,&BETA,GX_t,&d);

    // Rescale all eigenvalues by scaling factor
    F77_NAME(drscl)(&d2,&scale_factor,GX_t,&INC1);
    F77_NAME(drscl)(&d2,&scale_factor,GS_t,&INC1);
}

// X is untouched
// REQUIRES dwork 2d^2 + 3d + 1 (assuming dsyevd only needs 2d+1)
// REQUIRES iwork 1.
void smoothed_objective_nok(problem_instance* prob, double* X, double* lambda_min_tp1,
                            double* obj_tp1, workspace* work){
    // Grab Problem Parameters
    int d = prob -> d;
    double mu_n = prob -> mu_n;

    // Local Vars
    double* GS_t;
    double* d2_tmp;
    double S_min;
    double X_min;
    double lambda_min_n;
    double lambda_min;
    double obj_value;
    double* X_eigs;
    int X_eigs_idx;
    int d2 = d*d;
    int lapack_info;
    
    X_eigs_idx = work -> dsyevd_ldwork_N;
    X_eigs = (work -> dwork) + X_eigs_idx;
    GS_t = X_eigs + d;
    d2_tmp = GS_t + d2;

    // GS_t is just current value of X
    memcpy(GS_t,X,d2*sizeof(double));
    // eigenvectors will overwrite input matrix, so GX_t becomes V, eigenvectors
    memcpy(d2_tmp,X,d2*sizeof(double));
    F77_CALL(dsyevd)(&JOBZN,&UPLO,&d,d2_tmp,&d,X_eigs,work->dwork,&(work->dsyevd_ldwork_N),
                work->iwork,&(work->dsyevd_liwork_N),&lapack_info);

    X_min = min_array(d,X_eigs);
    S_min = min_array(d2,GS_t);
    lambda_min = S_min < X_min ? S_min : X_min;
    lambda_min_n = -1*lambda_min;
    *lambda_min_tp1 = lambda_min;


    // subtract lambda_min from all eigenvalues
    daps(X_eigs,1,lambda_min_n,d);
    daps(GS_t,1,lambda_min_n,d2);

    // Divide both eigenvalues by mu
    F77_NAME(drscl)(&d,&mu_n,X_eigs,&INC1);
    F77_NAME(drscl)(&d2,&mu_n,GS_t,&INC1);

    // Exponentiate eigenvalues
    dvexp(X_eigs,d);
    dvexp(GS_t,d2);

    // Compute Objective Value
    obj_value = dsumv(X_eigs,d);
    obj_value = obj_value + dsumv(GS_t,d2);
    obj_value = mu_n*log(obj_value) + lambda_min;

    *obj_tp1 = obj_value;
}
