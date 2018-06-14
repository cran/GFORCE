#include "float.h"
#include "math.h"
#include "time.h"
#include "R.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "util_FORCE.h"
#include "util_mops.h"

//prototypes
void project_C_perpendicular(problem_instance* prob, double* GX_t, double* GS_t, workspace* work);
void project_C_perpendicular_nok(problem_instance* prob, double* GX_t, double* GS_t, workspace* work);

// REQUIRES work -> dwork be of length at least d^2 + 7d + 6?
// REQUIRES work -> iwork be of length at least d+2?
// b requires d+2
// M requires d^2 + 4d + 4
// S_rsums requires d
// X_rsums requires d
// Return value is stored as X_t
void C_perp_update(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace* work) {
    double g_norm;
    int d = prob -> d;
    int d2 = d*d;

    // compute projection -- gradient is partially overwritten
    project_C_perpendicular(prob, GX_t, GS_t, work);

    // compute norm and rescale
    g_norm = F77_NAME(dnrm2)(&d2,GX_t,&INC1);
    g_norm = GNORM_FACTOR*g_norm; // need to account for slack variables
    F77_NAME(drscl)(&d2,&g_norm,GX_t,&INC1);

    // update and return
    F77_NAME(daxpy)(&d2,&alpha,GX_t,&INC1,X_t,&INC1);
}

// Return value is stored in Z_proj
// Z IS NOT ALTERED
void project_E(problem_instance* prob, double* Z, double lmin, double* Z_proj){
    int d = prob -> d;
    int d2 = d*d;
    double* E = prob -> E;
    double scale_factor = 1/(1-lmin);

    memcpy(Z_proj,Z,d2*sizeof(double));
    F77_NAME(dscal)(&d2,&scale_factor,Z_proj,&INC1);
    scale_factor = 1 - scale_factor;
    F77_NAME(daxpy)(&d2,&scale_factor,E,&INC1,Z_proj,&INC1);
}


// OUTPUT STORED IN GX_T
void project_C_perpendicular(problem_instance* prob, double* GX_t, double* GS_t, workspace* work){
    // Local Declarations
    int lapack_info;
    double scale_factor;
    double dtmp1, dtmp2;
    double* ptmp1;

    // Problem Instance Extraction, Workspace setup
    double* D_rsums = prob -> D_rsums;
    double* D = prob -> D;
    double TrD = prob -> TrD;
    double DTD = prob -> DTD;
    int d = prob -> d;
    int* ipiv = work -> iwork;
    int dp2 = d+2;
    int dp1 = d+1;
    int dp3 = d+3;
    int dp22 = dp2*dp2;
    int d2 = d*d;

    // Setup linear system, precompute
    double TrX = dtrace(GX_t,d);
    double TrS = dtrace(GS_t,d);
    double DTX = F77_NAME(ddot)(&d2,D,&INC1,GX_t,&INC1);
    double DTS = F77_NAME(ddot)(&d2,D,&INC1,GS_t,&INC1);
    double* M = work -> dwork;
    double* b = M + dp22;
    double* S_rsums = b + dp2;
    double* X_rsums = S_rsums + d;
    dcsum(GX_t,d,X_rsums);
    dcsum(GS_t,d,S_rsums);

    // M all ones
    ptmp1 = M;
    for(int i=0; i < dp22; i++){
        *ptmp1 = 1;
        ptmp1++;
    }
    //M <- M + p*diag(p+2)
    ptmp1 = M;
    for(int i=0; i < d; i++){
        *ptmp1 = dp1;
        ptmp1 = ptmp1 + dp3;
    }

    // M[1:p,p+2] <- rowSums(C)
    memcpy(M + dp1*dp2,D_rsums,d*sizeof(double));
    //b[1:p] <- rowSums(Z) + rowSums(S)
    dxpyez(d,X_rsums,S_rsums,b);
    //M[p+1,1:p] <- 2*rep(1,p)
    ptmp1 = M + d;
    for(int i=0; i < d; i++){
        *ptmp1 = 2;
        ptmp1 = ptmp1 + dp2;
    }
    // M[p+1,p+1] <- p
    *(M + d + d*dp2) = d;
    // M[p+1,p+2] <- TC
    *(M + d + (dp2-1)*dp2) = TrD;
    // b[p+1] <- sum(diag(Z)) + k3*sum(diag(S))
    *(b+d) = TrX + TrS;
    // M[p+2,1:p] <- 2*colSums(C)
    scale_factor = 2.0;
    ptmp1 = M+d+1;
    for(int i=0; i < d; i++){
        dtmp1 = (*D_rsums)*scale_factor;
        *ptmp1 = dtmp1;
        ptmp1 = ptmp1+dp2;
        D_rsums++;
    }
    // M[p+2,p+1] <- TC
    *(M + dp1 + d*dp2) = TrD;
    // M[p+2,p+2] <- sum(C*C)
    *(M + dp22 - 1) = DTD;
    // b[p+2] <- sum(C*Z) + k3*sum(C*S)
    *(b+dp1) = DTX + DTS;

    // b will store x, the solution
    F77_CALL(dgesv)(&dp2,&NRHS1,M,&dp2,ipiv,b,&dp2,&lapack_info);

    // multiply GX by 1/2 and store in GX
    scale_factor = 0.5;
    F77_NAME(dscal)(&d2,&scale_factor,GX_t,&INC1);

    // 0.5GX + 0.5GS -> GX
    F77_NAME(daxpy)(&d2,&scale_factor,GS_t,&INC1,GX_t,&INC1);

    // GX - 0.5*x[d+2]*D
    scale_factor = -0.5;
    dtmp1 = scale_factor * (*(b+d+1));
    F77_NAME(daxpy)(&d2,&dtmp1,D,&INC1,GX_t,&INC1);

    //GX - 0.5*x[d+1]*I
    dtmp1 = scale_factor * (*(b+d));
    ptmp1 = GX_t;
    for(int i=0; i < d; i++){
        dtmp2 = (*ptmp1) + dtmp1;
        *ptmp1 = dtmp2;
        ptmp1 = ptmp1 + dp1;
    }

    //GX -0.5*x[a]*R_a
    for(int a=0; a < d; a++){
        dtmp1 = scale_factor*(*(b+a));
        daps(GX_t + a*d,1,dtmp1,d);
        daps(GX_t + a,d,dtmp1,d);
    }
}


// REQUIRES work -> dwork be of length at least d^2 + 7d + 6?
// REQUIRES work -> iwork be of length at least d+2?
// b requires d+2
// M requires d^2 + 4d + 4
// S_rsums requires d
// X_rsums requires d
// Return value is stored as X_t
void C_perp_update_nok(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace* work) {
    double g_norm;
    int d = prob -> d;
    int d2 = d*d;

    // compute projection -- gradient is partially overwritten
    project_C_perpendicular_nok(prob, GX_t, GS_t, work);

    // compute norm and rescale
    g_norm = F77_NAME(dnrm2)(&d2,GX_t,&INC1);
    g_norm = GNORM_FACTOR*g_norm; // need to account for slack variables
    F77_NAME(drscl)(&d2,&g_norm,GX_t,&INC1);

    // update and return
    F77_NAME(daxpy)(&d2,&alpha,GX_t,&INC1,X_t,&INC1);
}


// OUTPUT STORED IN GX_T
void project_C_perpendicular_nok(problem_instance* prob, double* GX_t, double* GS_t, workspace* work){
    // Local Declarations
    int lapack_info;
    double scale_factor;
    double dtmp1;
    double* ptmp1;

    // Problem Instance Extraction, Workspace setup
    double* D_rsums = prob -> D_rsums;
    double* D = prob -> D;
    double DTD = prob -> DTD;
    int d = prob -> d;
    int* ipiv = work -> iwork;
    int dp2 = d+2;
    int dp1 = d+1;
    int dp12 = dp1*dp1;
    int d2 = d*d;

    // Setup linear system, precompute
    double DTX = F77_NAME(ddot)(&d2,D,&INC1,GX_t,&INC1);
    double DTS = F77_NAME(ddot)(&d2,D,&INC1,GS_t,&INC1);
    double* M = work -> dwork;
    double* b = M + dp12;
    double* S_rsums = b + dp1;
    double* X_rsums = S_rsums + d;
    dcsum(GX_t,d,X_rsums);
    dcsum(GS_t,d,S_rsums);

    // M all ones
    ptmp1 = M;
    for(int i=0; i < dp12; i++){
        *ptmp1 = 1;
        ptmp1++;
    }
    //M <- M + d*diag(d+1)
    ptmp1 = M;
    for(int i=0; i < d; i++){
        *ptmp1 = dp1;
        ptmp1 = ptmp1 + dp2;
    }

    // M[1:p,p+1] <- rowSums(D)
    memcpy(M + d*dp1,D_rsums,d*sizeof(double));
    //b[1:p] <- rowSums(Z) + rowSums(S)
    dxpyez(d,X_rsums,S_rsums,b);
    // M[p+1,1:p] <- 2*colSums(D)
    scale_factor = 2.0;
    ptmp1 = M+d;
    for(int i=0; i < d; i++){
        dtmp1 = (*D_rsums)*scale_factor;
        *ptmp1 = dtmp1;
        ptmp1 = ptmp1+dp1;
        D_rsums++;
    }
    // M[p+1,p+1] <- sum(D*D)
    *(M + dp12 - 1) = DTD;
    // b[p+1] <- sum(D*Z) + sum(D*S)
    *(b+d) = DTX + DTS;

    // b will store x, the solution
    F77_CALL(dgesv)(&dp1,&NRHS1,M,&dp1,ipiv,b,&dp1,&lapack_info);

    // multiply GX by 1/2 and store in GX
    scale_factor = 0.5;
    F77_NAME(dscal)(&d2,&scale_factor,GX_t,&INC1);

    // 0.5GX + 0.5GS -> GX
    F77_NAME(daxpy)(&d2,&scale_factor,GS_t,&INC1,GX_t,&INC1);

    // GX - 0.5*x[d+2]*D
    scale_factor = -0.5;
    dtmp1 = scale_factor * (*(b+d));
    F77_NAME(daxpy)(&d2,&dtmp1,D,&INC1,GX_t,&INC1);

    //GX -0.5*x[a]*R_a
    for(int a=0; a < d; a++){
        dtmp1 = scale_factor*(*(b+a));
        daps(GX_t + a*d,1,dtmp1,d);
        daps(GX_t + a,d,dtmp1,d);
    }
}
