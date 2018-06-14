#include "R.h"
#include "math.h"
#include "omp.h"
#include "FORCE.h"

void v_measure(double* restrict IPS,double* restrict n_xc_xd, int *dimension, double* restrict vm)
{
    int p = *dimension;
    /* return UPPER diagonal matrix of V measures b/c in COLMAJOR ORDER */
    for(int b=0; b < p; b++)
    {
        for(int a=0; a < b; a++)
        {
            double vmax = 0;
            for(int d=0; d < p; d++)
            {
                for(int c = 0; c < d; c++)
                {
                    if(c!=a && c != b && d!=a && d!=b && n_xc_xd[d*p+ c] != 0)
                    {
                        double vcd = (IPS[a*p+ c] + IPS[b*p+ d] - IPS[b*p+ c] - IPS[a*p+ d]) / n_xc_xd[d*p+ c];
                        double vcd_abs = (vcd) > 0 ? vcd : -1*vcd;
                        vmax = (vcd_abs) > (vmax) ? vcd_abs : vmax;
                    }
                }
            }
            vm[b*p+ a] = vmax;
        }
    }
}


/* return UPPER diagonal matrix of V measures b/c in COLMAJOR ORDER */
/* a is row and b is column */
void v_measure_par(double* restrict IPS,double* restrict n_xc_xd, int* dimension, double* restrict vm)
{
    int* tmp_ptr;
    int p = *dimension;
    int num_entries = p*(p-1) / 2;

    /* Need to iterate over only the upper triangle. That is indexing is like */
    /* for packed storage. A[i + j(j+1)/2 ], i == row, j== column */
    /* collapse nested loop, map indices */
    int* idx_map = (int*) R_alloc(2*num_entries,sizeof(int));
    tmp_ptr = idx_map;
    for(int b=0; b < p; b++){
        for(int a=0; a < b; a++){
            *tmp_ptr++ = a;
            *tmp_ptr++ = b;
        }
    }

    #pragma omp parallel for
    for(int entry_number=0; entry_number < num_entries; entry_number++)
    {
        double vmax = 0;
        int a = idx_map[2*entry_number];
        int b = idx_map[2*entry_number + 1];
        for(int d=0; d < p; d++)
        {
            for(int c = 0; c < d; c++)
            {
                if(c!=a && c != b && d!=a && d!=b && n_xc_xd[d*p+ c] != 0)
                {
                    double vcd = (IPS[a*p+ c] + IPS[b*p+ d] - IPS[b*p+ c] - IPS[a*p+ d]) / n_xc_xd[d*p+ c];
                    double vcd_abs = (vcd) > 0 ? vcd : -1*vcd;
                    vmax = (vcd_abs) > (vmax) ? vcd_abs : vmax;
                }
            }
        }
        vm[b*p+ a] = vmax;
    }
}
