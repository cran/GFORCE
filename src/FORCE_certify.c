#include "math.h"
#include "string.h"
#include "R.h"
#include "R_ext/Lapack.h"
#include "FORCE.h"
#include "util_FORCE.h"
#include "util.h"
#include "util_mops.h"

// FUNCTION PROTOTPES
void precompute_values(double* restrict Y_a_base, int* restrict group_sizes, double* restrict primal_value, double* restrict D, int* restrict ga_hat,
                        int d, int K, double* restrict dwork, int* restrict iwork);
void compute_Y_a(double* restrict Y_a_base, double* restrict Y_a_new, double Y_T, int* restrict ga_hat, int* restrict group_sizes, int d);

//KMEANS DUAL SOLUTION MIN PRIMAL FORM
//R ACCESS POINTS
void kmeans_dual_solution_primal_min_R(int* ga_hat, double* D, int* K_0, int *dimension, 
                                        double* eps1_0, double* eps2_0, double* Y_T_min_0, 
                                        double* Y_a_r, double* Y_T_r, int* feasible_r) {
    int d = *dimension;
    int K = *K_0;
    double Y_T_min = *Y_T_min_0;
    double eps1 = *eps1_0;
    double eps2 = *eps2_0;
    kmeans_dual_solution_primal_min(ga_hat,D,K,d,eps1,eps2,Y_T_min,Y_a_r,Y_T_r,feasible_r);
}

void kmeans_dual_solution_primal_min_nok_R(int* ga_hat, double* D, int *K_hat0, int *dimension, 
                                        double* eps1_0, double* Y_a_r, int* feasible_r) {
    int d = *dimension;
    int K_hat = *K_hat0;
    double eps1 = *eps1_0;
    kmeans_dual_solution_primal_min_nok(ga_hat,D,K_hat,d,eps1,Y_a_r,feasible_r);
}

//External C Access Points
void kmeans_dual_solution_primal_min(int* restrict ga_hat, double* restrict D, int K, int d, double eps1,
                                        double eps2, double Y_T_min, double* restrict Y_a_r,
                                        double* restrict Y_T_r, int* restrict feasible_r) {
    workspace work;
    problem_instance prob;
    prob.D = D;
    prob.K = K;
    prob.d = d;
    work.dwork = (double *) R_alloc((d*(d-1))/2 + 7*d -2,sizeof(double));
    work.iwork = (int *) R_alloc(d + 3*K + 3,sizeof(int));
    kmeans_dual_solution_impl(ga_hat,&prob,eps1,eps2,Y_T_min,Y_a_r,Y_T_r,feasible_r,&work);
}

void kmeans_dual_solution_primal_min_nok(int* restrict ga_hat, double* restrict D, int K_hat, int d,
                                        double eps1, double* restrict Y_a_r, int* restrict feasible_r) {
    workspace work;
    problem_instance prob;
    prob.D = D;
    prob.K = K_hat;
    prob.d = d;
    work.dwork = (double *) R_alloc(d*d + 7*d,sizeof(double));
    work.iwork = (int *) R_alloc(d + 3*K_hat + 3,sizeof(int));
    kmeans_dual_solution_nok_impl(ga_hat,&prob,eps1,Y_a_r,feasible_r,&work);
}

//Internal Access Point
//KMEANS DUAL SOLUTION MIN PRIMAL FORM
//ga_hat should be in "Standard Form" - group names are use 1...K
//REQUIRES d(d-1)/2 + 7d -2 length dwork
//REQUIRES d+3K+3 length iwork
//REQUIRES d length Y_a_r
void kmeans_dual_solution_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1, double eps2, 
                                double Y_T_min, double* restrict Y_a_r, double* restrict Y_T_r, int* restrict feasible_r,
                                workspace* restrict work) {
    // Local Variable Declarations
    double* Y_a_new;
    double* Y_a_base;
    double* Y_a_best;
    double* R;
    double* T_d; // for reduction to tridiagonal form
    double* T_e; // for reduction to tridiagonal form
    double* T_tau; // for reduction to tridiagonal form
    int* group_sizes;
    double dtmp1,Y_T_max,Y_T_best,em_min;
    int tmp1,tmp2;
    int iter_feasible,done,same_group; //boolean values
    double* D = prob -> D;
    int d = prob -> d;
    int K = prob -> K;

    // FOR LAPACK CALLS
    int lapack_info = 0;
    char lapack_format = 'U';

    // Initialization + Declaration

    int feasible = 0;
    double Y_T_new = 0;
    double primal_value = 0;
    eps1 = 1 - eps1;
    eps2 = -1*eps2;

    // Get memory from workspace for Y_a_base and group_sizes
    Y_a_base = work -> dwork;
    group_sizes = work -> iwork;

    // Pre-compute
    precompute_values(Y_a_base,group_sizes,&primal_value,D,ga_hat,d,K,Y_a_base+d,group_sizes+K+1);

    // Get rest of memory for working
    Y_a_new = Y_a_base + d;
    Y_a_best = Y_a_new + d;
    R = Y_a_best + d;
    T_d = R + d + (d*(d-1))/2;
    T_e = T_d + d;
    T_tau = T_e + d-1;

    //Binary Search For Y_T
    done = 0;
    Y_T_max = abs(primal_value);
    Y_T_best = Y_T_max;
    while(!done){
        //Update dual variables
        Y_T_new = (Y_T_max + Y_T_min)/2;
        compute_Y_a(Y_a_base, Y_a_new, Y_T_new, ga_hat, group_sizes, d);
        iter_feasible = 1;
        
        //iterate over off diagonal elements
        // This is will be stored in packed form
        // A[i + j(j+1)/2 ] i is row, j is column
        // iterate over column (b) then over row (a)
        // tmp2 stores index into R
        tmp2 = 0;
        for(int b = 0; b < d && iter_feasible; b++){
            for(int a = 0; a <= b && iter_feasible; a++) {
                tmp1 = b*d+a; // index into D
                dtmp1 = D[tmp1] + Y_a_new[a] + Y_a_new[b];
                same_group = ga_hat[a] == ga_hat[b];
                //checks dual feasibility of Y_ab for this choice of Y_T
                if(!same_group){
                    if(dtmp1 < 0.0){
                        iter_feasible = 0;
                    }
                    R[tmp2] = 0;
                } else {
                    R[tmp2] = dtmp1;
                }
                tmp2++;
            }
        }

        if(iter_feasible){
            //iterate over diagonal elements, and add Y_T x I
            for(int a = 0; a < d; a++){
                tmp1 = a + (a*(a+1))/2;
                dtmp1 = R[tmp1] + Y_T_new;
                R[tmp1] = dtmp1;
            }

            //find minimum eigenvalue -- first do transformation to tri-diagonal
            // R == ap, after call R has similarity transformation information
            F77_CALL(dsptrd)(&lapack_format,&d,R,T_d,T_e,T_tau,&lapack_info);

            //Get tridiagonal eigenvalues
            F77_CALL(dsterf)(&d,T_d,T_e,&lapack_info);

            em_min = min_array(d,T_d);

            //update bounds
            if(em_min > eps2){
                feasible = 1;
                Y_T_best = Y_T_new;
                Y_a_best = Y_a_new;
                Y_T_max = Y_T_new;
            } else{
                Y_T_min = Y_T_new;
            }

        } else{
            //update bounds
            Y_T_max = Y_T_new;
        }

        // Check for termination
        if(Y_T_min/Y_T_max > eps1){
            done = 1;
        }

    }

    // SET RETURN VALUES
    *feasible_r = feasible;
    *Y_T_r = Y_T_best;
    memcpy(Y_a_r,Y_a_best,d*sizeof(double));
}



//Internal Access Point
//KMEANS DUAL SOLUTION MIN PRIMAL FORM (K IS UNKNOWN)
//ga_hat should be in "Standard Form" - group names are use 1...K
//REQUIRES d(d-1)/2 + 7d -2 length dwork
//REQUIRES d+3K+3 length iwork
//REQUIRES d length Y_a_r
void kmeans_dual_solution_nok_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1,
                                   double* restrict Y_a_r, int* restrict feasible_r, workspace* restrict work) {
    // Local Variable Declarations
    double* R;
    double* T_d; // for reduction to tridiagonal form
    double* T_e; // for reduction to tridiagonal form
    double* T_tau; // for reduction to tridiagonal form
    int* group_sizes;
    double dtmp1,em_min;
    int tmp1;
    int same_group; //boolean values
    double* D = prob -> D;
    int d = prob -> d;
    int K_hat = prob -> K;

    // FOR LAPACK CALLS
    int lapack_info = 0;
    char lapack_format = 'U';

    // Initialization + Declaration
    int feasible = 1;
    double primal_value = 0;
    eps1 = -1*eps1;
    R = work -> dwork;
    T_d = R + d + (d*(d-1))/2;
    T_e = T_d + d;
    T_tau = T_e + d-1;

    // Get memory from workspace for group_sizes
    group_sizes = work -> iwork;

    // Compute Y_a (Y_T is always 0 when K is unknown)
    precompute_values(Y_a_r,group_sizes,&primal_value,D,ga_hat,d,K_hat,R,group_sizes+K_hat+1);

        
    //iterate over off diagonal elements
    // This is will be stored in packed form
    // A[i + j(j+1)/2 ] i is row, j is column
    // iterate over column (b) then over row (a)
    // tmp1 stores index into R
    tmp1 = 0;
    for(int b = 0; b < d && feasible; b++){
        for(int a = 0; a <= b && feasible; a++) {
            tmp1 = b*d+a; // index into D
            dtmp1 = D[tmp1] + Y_a_r[a] + Y_a_r[b];
            same_group = ga_hat[a] == ga_hat[b];
            //checks dual feasibility of Y_ab
            if(!same_group){
                if(dtmp1 < 0.0){
                    feasible = 0;
                }
                R[tmp1] = 0;
            } else {
                R[tmp1] = dtmp1;
            }
            tmp1++;
        }
    }

    if(feasible){
        //find minimum eigenvalue -- first do transformation to tri-diagonal
        // R == ap, after call R has similarity transformation information
        F77_CALL(dsptrd)(&lapack_format,&d,R,T_d,T_e,T_tau,&lapack_info);

        //Get tridiagonal eigenvalues
        F77_CALL(dsterf)(&d,T_d,T_e,&lapack_info);

        em_min = min_array(d,T_d);

        //check if psd
        if(em_min < eps1){
            feasible = 0;
        }

    }

    // SET RETURN VALUES
    *feasible_r = feasible;
}




// Computes Y_a_base, group_sums and group_sizes. Assumes these have already been allocated
// to the proper size. Also computes optval
//REQUIRES K+1 length dwork
//REQUIRES d+2K+2 length iwork
void precompute_values(double* restrict Y_a_base, int* restrict group_sizes, double* restrict primal_value,
                        double* restrict D, int* restrict ga_hat, int d, int K, double* restrict dwork, int* restrict iwork) {
    // Local Vars
    int tmp1,tmp2,tmp3,tmp4;
    double dtmp1;
    double* group_sums;
    int* group_tailp1_idx;
    int* group_start_idx; 
    int* group_idxs; // stores by group all idxs in that group
    double primal_value_0 = 0;

    // Allocate Memory
    group_tailp1_idx = iwork;
    group_start_idx = iwork+K+1;
    group_idxs = group_start_idx+K+1;
    group_sums = dwork;

    // Zero out
    for(int i=0; i < K+1; i++){
        group_sizes[i] = 0;
        group_sums[i] = 0;
    }

    // Group Sizes, zero out Y_a_base
    for(int i=0; i < d; i++){
        tmp1 = ga_hat[i];
        group_sizes[tmp1] = group_sizes[tmp1] + 1;
        Y_a_base[i] = 0;
    }

    // Group Bounds Initialization -- points to current last group member, first group member
    tmp1 = 0;
    for(int i=0; i < K+1; i++){
        group_tailp1_idx[i] = tmp1;
        group_start_idx[i] = tmp1;
        tmp1 += group_sizes[i];
    }

    // Group Membership
    for(int i=0; i < d; i++){
        tmp1 = ga_hat[i]; //store group membership of i
        tmp2 = group_tailp1_idx[tmp1];
        group_tailp1_idx[tmp1] = group_tailp1_idx[tmp1] + 1; //store location to add this index, then increment
        group_idxs[tmp2] = i;
    }

    // Row Sums -- can store temporarily in Y_a_base
    // Group Sum -- Can be done simultaneously
    for(int i=0; i < d; i++){
        tmp1 = ga_hat[i]; // group membership of i
        tmp2 = group_start_idx[tmp1]; //start idx of group
        tmp3 = group_tailp1_idx[tmp1]; // end idx
        dtmp1 = 0;
        for(int j=tmp2; j < tmp3; j++){
            //formula Y_a_base = -\frac{1}{|\Gs{i}|}\sum_{b \in \Gs{i}}D_{a,b}
            tmp4 = (group_idxs[j])*d + i; //Access is in column major form
            dtmp1 = dtmp1 + D[tmp4];
        }
        Y_a_base[i] = dtmp1;
        group_sums[tmp1] += dtmp1;
    }

    // Y_a_base
    // Primal Value -- can be computed at the same time
    // Primal Value -- \sum_{i} Y_a_base[i] == d^*
    for(int i=0; i < d; i++){
        tmp1 = ga_hat[i]; //get group of var i
        tmp2 = group_sizes[tmp1];
        dtmp1 = (group_sums[tmp1]/(2*tmp2) - Y_a_base[i])/tmp2;
        Y_a_base[i] = dtmp1;
        primal_value_0 += dtmp1;
    }
    //need to mult primal value by 2
    primal_value_0 = 2*primal_value_0;

    //Set Return Values
    *primal_value = primal_value_0;
}


void compute_Y_a(double* restrict Y_a_base, double* restrict Y_a_new, double Y_T, int* restrict ga_hat, int* restrict group_sizes, int d) {
    int tmp1;
    for(int i=0; i < d; i++) {
        tmp1 = ga_hat[i];
        Y_a_new[i] = Y_a_base[i] - Y_T/(2*group_sizes[tmp1]);
    }
}
