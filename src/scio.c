#include "R.h"
#include "math.h"
#include "string.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"
#include "glatent_inference.h"

static const char JOBZV = 'V';
static const char JOBZN = 'N';
static const char UPLO = 'L';
static const char SIDE_L = 'L';
static const char SIDE_R = 'R';
static const char TRANS_N = 'N';
static const char TRANS_T = 'T';
static const double ALPHA = -1.0;
static const double BETA = 0.0;
static const int INC1 = 1;
static const int NRHS1 = 1;


void scio_column_R(double* restrict C, int* d0, int* k0,double* restrict theta_k,
                    double* lambda0, double* eps0, int* max_iter0) {
    int d = *d0;
    int k = *k0;
    double lambda = *lambda0;
    double eps = *eps0;
    int max_iter = *max_iter0;

    double* dwork = (double *) R_alloc(d,sizeof(double));
    int* active_idx = (int *) R_alloc(d,sizeof(int));
    scio_column(C,d,k,theta_k,lambda,eps,max_iter,active_idx,dwork);
}


// requires active_idx of length d
// requires dwork of size d
void scio_column(double* restrict C, int d, int k, double* restrict theta_k,
                double lambda, double eps, int max_iter,int* restrict active_idx, double* restrict dwork) {

    //local vars
    double eps_act,delta,delta_t,theta_t_i,theta_tp1_i,dtmp,dtmp2;
    int done_eps;
    int* pitmp;
    double* pdtmp;
    double* nCtheta;
    int t;
    int active_set_mode;

    eps_act = eps*50; //arbitrary constant. 50 is same choice as scio package...
    done_eps = 0;
    t = 1;
    active_set_mode = 0;

    // init theta_k
    pdtmp = theta_k;
    for(int i=0; i < d; i++) {
        *pdtmp = 0;
        pdtmp++;
    } 
    theta_k[k] = C[k*d + k];

    //init nCtheta
    //output is in pdtmp
    nCtheta = dwork;
    // memcpy(pdtmp,theta_k,d*sizeof(double));
    F77_NAME(dgemv)(&TRANS_N,&d,&d,&ALPHA,C,&d,theta_k,&INC1,&BETA,nCtheta,&INC1);

    // init active indices
    pitmp = active_idx;
    for(int i=0; i < d; i++) {
        *pitmp = 0;
        pitmp++;
    }

    // main loop, each outer iteration loops over all indices in theta_k
    while(t < max_iter && done_eps == 0) {
        delta_t = 0;
        //inner loop
        for(int i=0; i < d; i++) {
            //only update if active_set_mode is 0 or active set mode is 1 and
            //the index is selected as active
            if(active_set_mode == 0 || active_idx[i] == 1) {
                theta_t_i = theta_k[i];
                theta_tp1_i = nCtheta[i];
                dtmp = C[i*d + i];
                theta_tp1_i = theta_tp1_i + dtmp*theta_t_i;
                
                // need to add 1 if diagonal element
                if(i == k){
                    theta_tp1_i = theta_tp1_i + 1;
                }

                // soft threshold theta_tp1_i
                dtmp = fabs(theta_tp1_i) - lambda;
                dtmp2 = theta_tp1_i >= 0 ? 1 : -1;

                if(dtmp > 0) {
                    theta_tp1_i = dtmp2*dtmp;
                    dtmp = C[i*d + i];
                    theta_tp1_i = theta_tp1_i / dtmp; //soft threshold and divide
                    // make the update
                    theta_k[i] = theta_tp1_i;
                } else{
                    theta_k[i] = 0;
                }

                // if theta_k[i] changes, need to update deltas and nCtheta
                if(theta_t_i != theta_k[i]) {
                    delta = theta_k[i] - theta_t_i;
                    dtmp = fabs(delta);
                    delta_t = delta_t > dtmp ? delta_t : dtmp;
                    for(int j=0; j < d; j++) {
                        dtmp = nCtheta[j] - delta*(C[i*d + j]);
                        nCtheta[j] = dtmp;
                    }
                }
            }
        }

        //update active sets
        if(active_set_mode == 1) {
            if(delta_t < eps) active_set_mode = 0;
        } else {
            if(delta_t < eps) done_eps = 1;
            if(delta_t < eps_act) {
                active_set_mode = 1;
                pitmp = active_idx;
                pdtmp = theta_k;

                //select active indices i.e. nonzero idxs
                for(int j=0; j < d; j++){
                    if(*pdtmp != 0.0) {
                        *pitmp = 1;
                    } else {
                        *pitmp = 0;
                    }
                    pitmp++;
                    pdtmp++;
                }
            }
        }

        t++;
    }

}
