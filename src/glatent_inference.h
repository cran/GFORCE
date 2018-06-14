#ifndef GLATENT_INFERENCE_H
#define GLATENT_INFERENCE_H

// scio estimator
void scio_column(double* restrict C, int d, int k, double* restrict theta_k,
                double lambda, double eps, int max_iter,int* restrict active_idx, double* restrict dwork);

void scio_column_R(double* restrict C, int* d0, int* k0,double* restrict theta_k,
                    double* lambda0, double* eps0, int* max_iter0);

#endif
