#include "float.h"
#include "math.h"
#include "R.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "util_mops.h"


// VECTOR, MATRIX OPS
void daps(double* restrict A, int inc_A, double c, int d){
    double dtmp;
    if(inc_A == 1){
        for(int i=0; i < d; i++){
            A[i] = A[i] + c;
        }
    } else{
        for(int i=0; i < d; i++){
            dtmp = *A + c;
            *A = dtmp;
            A = A + inc_A;
        }
    }
}

// performs y = a*x + b*y
void daxpby(double a, double* restrict X, double b, double* restrict Y, int d){
    double dtmp1;
    for(int i=0; i < d; i++){
        dtmp1 = a * (*X) + b * (*Y);
        *Y = dtmp1;
        X++;
        Y++;
    }
}

void dsmtd(double* restrict A, double* restrict B, int d, const char side){
    int stride = 1;
    if(side == 'R'){
        for(int i=0; i <d; i++){
            F77_NAME(dscal)(&d,B+i,A + d*i,&stride);
        }
    }else {
        stride = d;
        for(int i=0; i <d; i++){
            F77_NAME(dscal)(&d,B+i,A + i,&stride);
        }
    }
}

void dvexp(double* restrict A, int d) {
    double tmp;
    for(int i = 0; i < d; i++){
        tmp = *A;
        *A = exp(tmp);
        A++;
    }
}

double dsumv(double* restrict A, int d){
    double acc = 0;
    for(int i=0; i < d; i++){
        acc = acc + *A;
        A++;
    }
    return acc;
}

double dtrace(double* restrict A, int d){
    double trace = 0;
    for(int i=0; i < d; i++){
        trace = trace + *A;
        A = A + d + 1;
    }
    return trace;
}


void dcsum(double* restrict A, int d, double* restrict A_csums){
    // zero A_csums
    double* tmp_ptr = A_csums;
    double dtmp;
    for(int i=0; i <d; i++){
        *tmp_ptr = 0;
        tmp_ptr++;
    }
    
    for(int i=0; i < d; i++){
        for(int j=0; j < d; j++){
            dtmp = *A_csums + *A;
            *A_csums = dtmp;
            A++;
        }
        A_csums++;
    }
}

void dxpyez(int d, double* restrict X, double* restrict Y, double* restrict Z){
    double dtmp1;
    for(int i=0; i < d; i++){
        dtmp1 = *X + *Y;
        *Z = dtmp1;
        X++;
        Y++;
        Z++;
    }
}

// Compute x = tr(A * B(G))
// Assumes clusters labeled 1..K or 0..K-1
// REQUIRES work -> dwork be of length K
// REQUIRES work -> iwork be of length d+3K+3
double dabgtp(double* restrict A, int* restrict G, int d, int K,int* iwork,double* dwork) {
    // Local Vars
    double opt_val = 0.0;
    int tmp1,tmp2,tmp3,tmp4;
    double dtmp1;
    int Kp1 = K+1;

    // Get Memory From Workspace
    double* group_sums = dwork;
    int* group_sizes = iwork;
    int* group_tailp1_idx = group_sizes + Kp1;
    int* group_start_idx = group_tailp1_idx + Kp1; 
    int* group_idxs = group_start_idx + Kp1; // stores by group all idxs in that group

    // Zero out
    for(int i=0; i < Kp1; i++){
        group_sizes[i] = 0;
        group_sums[i] = 0;
    }

    // Group Sizes, zero out Y_a_base
    for(int i=0; i < d; i++){
        tmp1 = G[i];
        group_sizes[tmp1] = group_sizes[tmp1] + 1;
    }

    // Group Bounds Initialization -- points to current last group member, first group member
    tmp1 = 0;
    for(int i=0; i < Kp1; i++){
        group_tailp1_idx[i] = tmp1;
        group_start_idx[i] = tmp1;
        tmp1 += group_sizes[i];
    }

    // Group Membership
    for(int i=0; i < d; i++){
        tmp1 = G[i]; //store group membership of i
        tmp2 = group_tailp1_idx[tmp1];
        group_tailp1_idx[tmp1] = tmp2 + 1; //store location to add this index, then increment
        group_idxs[tmp2] = i;
    }

    // Group Sums
    for(int i=0; i < d; i++){
        tmp1 = G[i]; // group membership of i
        tmp2 = group_start_idx[tmp1]; //start idx of group
        tmp3 = group_tailp1_idx[tmp1]; // end idx

        for(int j=tmp2; j < tmp3; j++){
            tmp4 = (group_idxs[j])*d + i; //Access is in column major form
            dtmp1 = group_sums[tmp1] + A[tmp4];
            group_sums[tmp1] = dtmp1;
        }
    }

    // Objective Value
    opt_val = 0;
    for(int i=0; i < Kp1; i++){
        tmp1 = group_sizes[i];
        if(tmp1 > 0){
            opt_val = opt_val + (group_sums[i] / ((double) tmp1));
        }
    }

    return opt_val;

}
