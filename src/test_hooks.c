#include "FORCE.h"
#include "R.h"
#include "R_ext/Lapack.h"
#include "test_hooks.h"
#include "util.h"
#include "util_mops.h"
#include "util_FORCE.h"


void test_smoothed_gradient_S_base(double* X, double* E,double* GS_t, int* d, double* S_min_r) {
    int d0 = *d;
    smoothed_gradient_S_base(X, E, GS_t, d0, S_min_r);
}


void test_smoothed_gradient_X_base(double* X, double* ESI, double* GX_t, int* d0, int* K0, double* eig_vals){
    workspace work;
    int d = *d0;
    int K = *K0;
    int X_eigs_idx;
    int tmp;
    int d2 = d*d;

    double* d2_tmp = (double *) R_alloc(d2,sizeof(double));
    allocate_workspace_FORCE(d,K,&work);

    X_eigs_idx = work.dsyevd_ldwork;
    tmp = 2*d2;
    X_eigs_idx = X_eigs_idx > tmp ? X_eigs_idx : tmp;

    smoothed_gradient_X_base(X, ESI, GX_t, d2_tmp, &work, d,X_eigs_idx);
    memcpy(eig_vals,work.dwork + X_eigs_idx,d*sizeof(double));
}

void test_smoothed_gradient(double* X, double* E,double* ESI,int* d0, int* K0,
                            double* mu0,double* GX_t,double* GS_t){
    workspace work;
    problem_instance prob;
    int d = *d0;
    int K = *K0;

    prob.E = E;
    prob.ESI = ESI;
    prob.d = d;
    prob.K = K;
    prob.mu_n = -1 * (*mu0);

    allocate_workspace_FORCE(d,K,&work);
    smoothed_gradient(&prob,X,GX_t,GS_t,&work);
}

void test_smoothed_objective(double* X, double* E,double* ESI,int* d0, int* K0,
                            double* mu0, double* lambda_min, double* obj_val){
    workspace work;
    problem_instance prob;
    int d = *d0;
    int K = *K0;

    prob.E = E;
    prob.ESI = ESI;
    prob.d = d;
    prob.K = K;
    prob.mu_n = -1 * (*mu0);

    allocate_workspace_FORCE(d,K,&work);
    smoothed_objective(&prob,X,lambda_min,obj_val,&work);
}

void test_project_C_perpendicular(double* D, int* d0, int* K0, double* GX_t, double* GS_t){
    workspace work;
    problem_instance prob;
    int d = *d0;
    int K = *K0;
    double mu = 1;
    double* E = NULL;
    double* ESI = NULL;

    allocate_workspace_FORCE(d,K,&work);
    initialize_problem_instance(D,E,ESI,mu,d,K,&prob);
    project_C_perpendicular(&prob,GX_t,GS_t,&work);
}

void test_project_E(double* E, double* Z, int* d, double* lmin, double* Z_proj){
    problem_instance prob;
    prob.d = *d;
    prob.E = E;
    project_E(&prob,Z,*lmin,Z_proj);

}

void test_clust_to_opt_val(double* D, int* d, int* K, int* clusters, double* opt_val_r){
    double opt_val;
    problem_instance prob;
    workspace work;
    prob.D = D;
    prob.d = *d;
    prob.K = *K;
    allocate_workspace_FORCE(*d,*K,&work);
    opt_val = clust_to_opt_val(&prob,clusters,&work);
    *opt_val_r = opt_val;
}

void test_project_C_perpendicular_nok(double* D, int* d0, double* GX_t, double* GS_t){
    workspace work;
    problem_instance prob;
    int d = *d0;
    double mu = 1;
    double* E = NULL;
    double* ESI = NULL;

    allocate_workspace_FORCE(d,d,&work);
    initialize_problem_instance(D,E,ESI,mu,d,d,&prob);
    project_C_perpendicular_nok(&prob,GX_t,GS_t,&work);
}


void test_daps(double* A,double* c, int* d){
    daps(A,1,*c,*d);
}

void test_dsmtd(double* A,double* B, int* d, char** side){
    dsmtd(A,B,*d,**side);
}

void test_dvexp(double* A, int* d){
    dvexp(A,*d);
}

void test_dsumv(double* A, int* d, double* sum_r){
    double tmp;
    tmp = dsumv(A,*d);
    *sum_r = tmp;
}

void test_dtrace(double* A, int* d, double* trace_r){
    double tmp;
    tmp = dtrace(A,*d);
    *trace_r = tmp;
}

void test_dcsum(double* A, int* d, double* A_csums){
    dcsum(A,*d,A_csums);
}

void test_dxpyez(int* d, double* X, double* Y, double* Z){
    dxpyez(*d,X,Y,Z);
}
