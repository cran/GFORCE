#include "R.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"
#include "FORCE.h"
#include "util.h"
#include "util_mops.h"


// void hclust(double* dist,int n,int m,int* agglomerate_idx_1, int* agglomerate_idx_2, double* agglomerate_dmin,double* dwork,int* iwork);

int L_curve_criterion(double* vals,int n);
void hclust_agg2clust(int* ag1,int* ag2,int n,int K,int* clusters);
void hclust_distances(double* X,int n, int m,double* dists,double* dwork);

void hclust_agglomerate_R(double* data, int* n0, int* agglomerate_idx_1, int* agglomerate_idx_2, double* agglomerate_dmin){
    hclust_agg_t hclust_sol;
    int n = *n0;
    int* iwork;
    double* dwork;
    hclust_sol.agg_idx_min = agglomerate_idx_1;
    hclust_sol.agg_idx_max = agglomerate_idx_2;
    hclust_sol.agg_dist = agglomerate_dmin;
    dwork = (double *) R_alloc(n,sizeof(double));
    iwork = (int *) R_alloc(2*n,sizeof(int));

    hclust_agglomerate(data,n,&hclust_sol,dwork,iwork);
}

void hclust_R(double* data, int* n0, int* clusters,int* K, double* MSE){
    hclust_t hclust_sol;
    int n = *n0;
    int* iwork;
    double* dwork;
    hclust_sol.MSE = MSE;
    hclust_sol.clusters = clusters;
    dwork = (double *) R_alloc(2*n,sizeof(double));
    iwork = (int *) R_alloc(7*n + 3,sizeof(int));

    hclust(data,n,&hclust_sol,dwork,iwork);

    // set return value for K
    *K = hclust_sol.K;
}

//requires dwork of length at least 3*d^2
//requires iwork of length at least 7d + 3
void hclust_FORCE(double* D,int d,hclust_t* hclust_sol,workspace* work) {
    double* dwork = work -> dwork;
    double* dists = dwork;
    dwork = dwork + d*d;

    hclust_distances(D,d,d,dists,dwork);

    hclust_sol -> MSE = dwork;
    dwork = dwork + d;
    int* iwork = work -> iwork;

    hclust(dists,d,hclust_sol,dwork,iwork);
}


// returns hclust_t
//requires dwork of length at least 2n
//requires iwork of length at least 7n + 3
void hclust(double* dists,int n,hclust_t* hclust_sol,double* dwork,int* iwork) {
    hclust_agg_t hclust_agg;
    int* clusters;
    double* mse;
    int* ag1;
    int* ag2;
    int idx1,idx2,clust1,clust2;

    // Step 0 - Init
    mse = hclust_sol -> MSE;
    hclust_agg.agg_dist = dwork;
    dwork = dwork + n;
    hclust_agg.agg_idx_min = iwork;
    iwork = iwork + n;
    hclust_agg.agg_idx_max = iwork;
    iwork = iwork + n;
    hclust_agg.n = n;
    clusters = iwork;
    iwork = iwork + n;

    for(int i=0; i < n; i++){
        clusters[i] = i;
    }

    // Step 1 - Perform Agglomeration
    hclust_agglomerate(dists,n,&hclust_agg,dwork,iwork);
    ag1 = hclust_agg.agg_idx_min;
    ag2 = hclust_agg.agg_idx_max;

    // Step 2 - Compute all MSEs
    int k = n;
    while(k > 1) {
        mse[k-1] = dabgtp(dists,clusters,n,k,iwork,dwork);

        // Update clusterings based on agglomeration
        idx1 = ag1[n-k];
        idx2 = ag2[n-k];
        clust1 = clusters[idx1];
        clust2 = clusters[idx2];

        // Relabel clusters
        // name clust2 with label clust1
        for(int i=0; i <n; i++){
            if(clusters[i] == clust2) {
                clusters[i] = clust1;
            }
        }

        // reorder cluster labels
        for(int i=0; i <n; i++){
            if(clusters[i] == k-1) {
                clusters[i] = clust2;
            }
        }

        // Decrement Cluster Number
        k--;
    }
    // calculate vals for last clustering
    mse[0] = dabgtp(dists,clusters,n,1,iwork,dwork);

    // Step 3 - Find minimum MSE and Estimate K
    k = L_curve_criterion(mse+1,n-1);
    k  = k + 2;
    hclust_sol -> K = k;

    // Step 4 - Get Clustering
    hclust_agg2clust(ag1,ag2,n,k,hclust_sol->clusters);
}


void hclust_agg2clust(int* ag1,int* ag2,int n,int K,int* clusters) {
    int idx1,idx2,clust1,clust2;

    for(int i=0; i < n; i++){
        clusters[i] = i;
    }


    for(int i=0; i < n-K; i++) {
        // Update clusterings based on agglomeration
        idx1 = ag1[i];
        idx2 = ag2[i];
        clust1 = clusters[idx1];
        clust2 = clusters[idx2];

        // Relabel clusters
        // name clust2 with label clust1
        for(int j=0; j <n; j++){
            if(clusters[j] == clust2) {
                clusters[j] = clust1;
            }
        }

        // reorder cluster labels
        for(int j=0; j <n; j++){
            if(clusters[j] == n-i-1) {
                clusters[j] = clust2;
            }
        }
    }
}

int L_curve_criterion(double* vals,int n){
    int max_idx;
    double dmax,a_1,a_2,b_1,b_2,n_1,n_2,nnorm,c_1,c_2,cproj_l,cperp_1,cperp_2,dnew;
    a_1 = 1;
    a_2 = vals[0];
    b_1 = n;
    b_2 = vals[n-1];
    n_1 = b_1 - a_1;
    n_2 = b_2 - a_2;
    nnorm = sqrt(pow(n_1,2) + pow(n_2,2));
    n_1 = n_1 / nnorm;
    n_2 = n_2 / nnorm;

    dmax = -DBL_MAX;
    max_idx = -1;
    for(int i=0; i < n; i++) {
        c_1 = i;
        c_2 = vals[i];
        cproj_l = (a_1 - c_1)*n_1 + (a_2 - c_2)*n_2;
        cperp_1 = (a_1 - c_1) - n_1*cproj_l;
        cperp_2 = (a_2 - c_2) - n_2*cproj_l;
        dnew = cperp_1*cperp_1 + cperp_2*cperp_2;

        if(dnew > dmax){
            dmax = dnew;
            max_idx = i;
        }
    }

    return(max_idx);
}



//requires dwork of length at least n
//requires iwork of length at least 2n
//WARNING: dist is overwritten
void hclust_agglomerate(double* dist,int n,hclust_agg_t* hclust_sol,double* dwork,int* iwork) {
    // Step 0 - Declarations
    int* nn_idx; //nearest neighbor idx
    int* active;
    int* agglomerate_idx_1;
    int* agglomerate_idx_2;
    double* agglomerate_dmin;
    double* nn_dist; //nearest neighbor distance
    int n2,itmp1,itmp2,idx_min,idx_max,dmin_idx1,dmin_idx2,num_clust;
    double dtmp1,dtmp2,dmin;

    // Step 1 - Init
    nn_dist = dwork;
    nn_idx = iwork;
    active = iwork + n;
    n2 = n*n;
    agglomerate_idx_1 = hclust_sol -> agg_idx_min;
    agglomerate_idx_2 = hclust_sol -> agg_idx_max;
    agglomerate_dmin = hclust_sol -> agg_dist;

    // Step 2 - Active Flags

    for(int i=0; i < n; i++){
        active[i] = 1;
    }


    // Step 3 - Initial nearest neighbors
    // This finds nearest neighbor with larger index
    for(int i=0; i < n-1;i++) {
        dtmp1 = DBL_MAX; //min distance
        itmp1 = -1; //min dist idx
        for(int j=i+1; j < n; j++) {
            itmp2 = i*n + j; //offset
            if(dtmp1 > dist[itmp2]) {
                dtmp1 = dist[itmp2];
                itmp1 = j;
            }

        }
        nn_idx[i] = itmp1;
        nn_dist[i] = dtmp1;
    }

    // Step 4 - 
    num_clust = n;

    while(num_clust > 1) {

        // find closest distance
        dmin = DBL_MAX;
        dmin_idx1 = -1;
        dmin_idx2 = -1;   
        for(int i=0; i < n-1; i++){
            if(active[i] == 1 && nn_dist[i] < dmin) {
                dmin = nn_dist[i];
                dmin_idx1 = i;
                dmin_idx2 = nn_idx[i];
            }
        }
        

        // add to list of agglomerations
        idx_min = dmin_idx1 > dmin_idx2 ? dmin_idx2 : dmin_idx1; // smaller idx
        idx_max = dmin_idx1 > dmin_idx2 ? dmin_idx1 : dmin_idx2; // larger idx
        
        agglomerate_idx_1[n-num_clust] = idx_min;
        agglomerate_idx_2[n-num_clust] = idx_max;
        agglomerate_dmin[n-num_clust] = dmin;

        // deactivate larger idx (bc we found nearest neighbors to the right)
        active[idx_max] = 0;

        // update distances from the new merged cluster
        dmin = DBL_MAX;
        for(int i=0; i < n; i++) {
            if(active[i] == 1 && i != idx_min) {
                //update dist matrix
                dtmp1 = dist[idx_min*n + i];
                dtmp2 = dist[idx_max*n + i];
                dtmp1 = dtmp1 > dtmp2 ? dtmp1 : dtmp2;
                dist[idx_min*n+ i] = dtmp1;
                dist[i*n+ idx_min] = dtmp1;

                //update dmin for idx_min to the right
                //check if need to update nn_idx,nn_dist for i
                dtmp1 = dist[idx_min*n + i];
                dtmp2 = nn_dist[i];
                if(idx_min < i && dtmp1 < dmin) {
                    dmin = dtmp1;
                    dmin_idx2 = i;
                } else if(idx_min > i && dtmp1 < dtmp2) {
                    nn_dist[i] = dtmp1;
                    nn_idx[i] = idx_min;
                }
            }
        }

        // update nearest neighbor for idx_min
        nn_idx[idx_min] = dmin_idx2;
        nn_dist[idx_min] = dmin;

        //update nearest neighbor to the right as required
        for(int i=0; i < n-1; i++){
            if(active[i] == 1 && (nn_idx[i] == idx_min || nn_idx[i] == idx_max)) {
                dmin = DBL_MAX;
                dmin_idx1 = -1;
                for(int j = i+1; j < n; j++) {
                    dtmp1 = dist[i*n + j];
                    if(active[j] == 1 && dtmp1 < dmin) {
                        dmin = dtmp1;
                        dmin_idx1 = j;
                    }
                }
                nn_idx[i] = dmin_idx1;
                nn_dist[i] = dmin;
            }
        }


        // decrement cluster number for next agglomeration
        num_clust = num_clust - 1;
    }

}


// construct distances matrix
// n x m matrix, n data points
// output is n x n 
// requires dwork = n^2 + n*m
void hclust_distances(double* X,int n, int m,double* dists,double* dwork) {
    double dtmp1;
    double* ips = dwork;
    dwork = dwork + n*n;
    double* X_copy = dwork;
    memcpy(X_copy,X,n*n*sizeof(double));

    //neither matrix is symmetric, so need to use dgemm
    F77_NAME(dgemm)(&TRANS_N,&TRANS_T,&n,&n,&m,&ALPHA,X,&n,X_copy,&n,&BETA,ips,&n);

    //fill dists matrix
    for(int i=0; i < n; i++) {
        for(int j=0; j < n; j++) {
            dtmp1 = ips[i*(n+1)]; //gets ith diagonal element
            dtmp1 = dtmp1 + ips[j*(n+1)]; //gets jth diagonal element
            dtmp1 = dtmp1 - 2*(ips[i*n + j]);
            dists[i*n+j] = sqrt(dtmp1);
        }
    }
}


