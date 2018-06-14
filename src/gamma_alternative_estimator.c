#include "R.h"
#include "math.h"
#include "omp.h"
#include "FORCE.h"

void gamma_alternative_estimator_R(double* restrict IPS,double* restrict ips_diag_sqrt, int* dimension, 
                                    double* scaling, int* restrict nes, double* restrict gamma_hat)
{
    int d = *dimension;
    double scale_factor = *scaling;
    double* work = (double *) R_alloc(d*d,sizeof(double));
    gamma_alternative_estimator(IPS,ips_diag_sqrt,d,scale_factor,nes,gamma_hat,work);
}

void gamma_alternative_estimator(double* restrict IPS,double* restrict ips_diag_sqrt, int d, double scaling, 
                                    int* restrict nes, double* restrict gamma_hat, double* restrict ne_meas)
{
    /* First find for all pairs a and b, max|<x_a-x_b,x_c/||x_c||>|. Use col-major order. */
    for(int b=0; b < d; b++)
    {
        for(int a=0; a < b; a++)
        {
            double max_c_val = -1;
            for(int c=0; c < d; c++)
            {
                if(c!=a && c != b)
                {
                    double new_c_val = (IPS[a*d+c]-IPS[b*d+c])/ips_diag_sqrt[c];
                    double c_val_abs = (new_c_val) > 0 ? new_c_val : -1*new_c_val;
                    max_c_val = (c_val_abs) > (max_c_val) ? c_val_abs : max_c_val;
                }
            }
            ne_meas[b*d+ a] = max_c_val;
            ne_meas[a*d+ b] = max_c_val;
        }
    }

    /* Get each variable's neighbor */
    for(int a=0; a < d; a++){
        int min_idx = -1;
        double min_val = 0;
        for(int b=0; b < d; b++){
            double new_val = ne_meas[a*d+b];
            if((a != b) && (new_val < min_val || min_idx == -1)){
                min_idx = b;
                min_val = new_val;
            }
        }
        nes[a] = min_idx;
        gamma_hat[a] = scaling*(IPS[a*d+a] - IPS[a*d+min_idx]);
    }
}


void gamma_alternative_estimator_par_R(double* restrict IPS,double* restrict ips_diag_sqrt, int* dimension,
                                        double* scaling, int* restrict nes, double* restrict gamma_hat)
{
    int d = *dimension;
    double scale_factor = *scaling;
    double* work = (double *) R_alloc(d*d,sizeof(double));
    int* idx_map = (int*) R_alloc(d*(d-1),sizeof(int));
    gamma_alternative_estimator_par(IPS,ips_diag_sqrt,d,scale_factor,nes,gamma_hat,work,idx_map);
}

/* par_idxs should have size d*(d-1) */
void gamma_alternative_estimator_par(double* restrict IPS,double* restrict ips_diag_sqrt, int d, double scaling,
                                        int* restrict nes, double* restrict gamma_hat, double* restrict ne_meas, int* restrict par_idxs)
{
    int* tmp_ptr;
    int num_entries = d*(d-1) / 2;
    /* need to populate par_idxs */
    tmp_ptr = par_idxs;
    for(int b=0; b < d; b++){
        for(int a=0; a < b; a++){
            *tmp_ptr++ = a;
            *tmp_ptr++ = b;
        }
    }

    /* First find for all pairs a and b, max|<x_a-x_b,x_c/||x_c||>|. Use par_idxs mapping. */
    #pragma omp parallel for
    for(int entry_number=0; entry_number < num_entries; entry_number++)
    {
        int a = par_idxs[2*entry_number];
        int b = par_idxs[2*entry_number + 1];
        double max_c_val = -1;
        for(int c=0; c < d; c++)
        {
            if(c!=a && c != b)
            {
                double new_c_val = (IPS[a*d+c]-IPS[b*d+c])/ips_diag_sqrt[c];
                double c_val_abs = (new_c_val) > 0 ? new_c_val : -1*new_c_val;
                max_c_val = (c_val_abs) > (max_c_val) ? c_val_abs : max_c_val;
            }
        }
        ne_meas[b*d+ a] = max_c_val;
        ne_meas[a*d+ b] = max_c_val;
    }

    /* Get each variable's neighbor */
    #pragma omp parallel for
    for(int a=0; a < d; a++){
        int min_idx = -1;
        double min_val = 0;
        for(int b=0; b < d; b++){
            double new_val = ne_meas[a*d+b];
            if((a != b) && (new_val < min_val || min_idx == -1)){
                min_idx = b;
                min_val = new_val;
            }
        }
        nes[a] = min_idx;
        gamma_hat[a] = scaling*(IPS[a*d+a] - IPS[a*d+min_idx]);
    }
}
