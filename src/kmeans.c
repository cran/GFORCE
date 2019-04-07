#include "math.h"
#include "string.h"
#include "time.h"
#include "R.h"
#include "omp.h"
#include "FORCE.h"
#include "util_FORCE.h"
#include "util.h"
#include "util_mops.h"

// Local Constant Strings
static const char* RANDOM_INIT = "random";
// static const char* RANDOM = "random";
static const char* KMEANS_PP_INIT = "kmeans++";
static const char* GIVEN_INIT = "kmeans++";

// Function prototypes
void lloyd_update_centers(double const* restrict D, double* restrict centers, int* restrict cluster_assignment,
                            int n, int m, int K, int* restrict iwork);
int lloyd_update_clusters(double const* restrict D, double* restrict centers, int* restrict cluster_assignment,
                            int n, int m, int K, double* restrict euclidean_distance_tmp);
int sample_discrete_distribution(double* restrict prob_dist,int n);
void update_min_distance(double const* restrict D, double* restrict min_center_distance, int new_center_idx,
                            int n, int m, double* restrict euclidean_distance_tmp);
void min_distance_to_probability(double* restrict min_distances, double* restrict prob_dist, int n);
double euclidean_distance(double const* restrict p1, double const* restrict p2, int m,double* restrict euclidean_distance_tmp);


// // Generic Kmeans -- specify initialization method
// void kmeans(double* D, int K, int n, int m, const char* init_mode,int* centers_init, int* cluster_assignment_r, double* centers_r){
//     if(strcmp(KMEANS_PP_INIT,init_mode)){
//         ;
//     } else if(strcmp(RANDOM_INIT,init_mode)){
//         ;
//     } else if(strcmp(GIVEN_INIT,init_mode)) {
//         ;
//     } else{
//         ;
//     }

// }


// R ACCESS POINT
void kmeans_pp_R(double* D, int* K0, int* n0, int* m0, int* cluster_assignment_r, double* centers_r, int* num_iters_R, double* time_R){
    int K = *K0;
    int n = *n0;
    int m = *m0;
    kmeans_pp(D,K,n,m,cluster_assignment_r,centers_r, num_iters_R, time_R);
}

// C ACCESS POINT
void kmeans_pp(double* D, int K, int n, int m, int* cluster_assignment_r, double* centers_r, int* num_iters_R, double* time_R) {
    workspace work;
    work.dwork = (double *) R_alloc(2*n + m*(K+1),sizeof(double));
    work.iwork = (int *) R_alloc(K+n,sizeof(int));
    kmeans_pp_impl(D,K,n,m,cluster_assignment_r,centers_r,num_iters_R,time_R,&work);
}


// void kmeans_impl(double* D, int K, int n, int m, int* cluster_assignment_r,
//                     double* centers_r, workspace* work) 


// INTERNAL ACCESS POINT
// Column major matrix layouts. Each column is a different data point.
// The matrix D should be m0 \times \n0 dimensional.
// NOT THREAD SAFE -- USES R INTERNALS RNG
// REQUIRES K+n length iwork
// REQUIRES 2n (prob dist,min_center_distance) + mK (centers) + m (euclidean_distance_tmp) = 2n + (K+1)m length dwork
void kmeans_pp_impl(double const* restrict D, int K, int n, int m, int* restrict cluster_assignment_r,
                    double* restrict centers_r,int* num_iters_R, double* time_R, workspace* work) {
    GetRNGstate();

    //////////////////////////////////////
    //// STEP 1 - K-means++ Initialization
    //////////////////////////////////////
    double* restrict prob_dist = work -> dwork;
    double* restrict min_center_distance = prob_dist + n;
    double* restrict centers = min_center_distance + n;
    double* restrict euclidean_distance_tmp = centers + m*K;
    int* restrict initial_centers = work -> iwork;
    int* restrict cluster_assignment = initial_centers + K;
    int tmp1;
    int num_iter = 0;
    double run_time = 0.0;
    struct timespec start_time, cur_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);


    //choose first center
    double q = unif_rand(); // random U[0,1]
    int c_idx = round((n+0.98-1)*q - 0.49); // approximately U[0,d-1]
    initial_centers[0] = c_idx;
    for(int i= 0; i < n; i++){
        if(i != c_idx){
            min_center_distance[i] = euclidean_distance(D+m*i,D+m*c_idx,m,euclidean_distance_tmp);
        } else {
            min_center_distance[i] = 0.0;
        }
    }
    min_distance_to_probability(min_center_distance,prob_dist,n);

    // generate centers 1..K-2
    for(int i=1; i < K-1; i++){
        tmp1 = sample_discrete_distribution(prob_dist,n);
        initial_centers[i] = tmp1;
        update_min_distance(D,min_center_distance,tmp1,n,m,euclidean_distance_tmp);
        min_distance_to_probability(min_center_distance,prob_dist,n);
    }
    // get last center (idx K-1)
    initial_centers[K-1] = sample_discrete_distribution(prob_dist,n);
    
    ///////////////////////////////////
    //// STEP 2 - K-means Clustering (Lloyd's Algorithm)
    ///////////////////////////////////
    int unchanged = 0;

    // initial clustering, zeros assignment to clusters, guarantees change in first iteration
    for(int k=0; k < K; k++){
        tmp1 = initial_centers[k];
        memcpy(centers + k*m,D+tmp1*m,m*sizeof(double));
    }
    for(int i=0; i < n; i++) cluster_assignment[i] = 0;

    // iterate
    while(!unchanged){
        // 1. Update Clusters
        unchanged = lloyd_update_clusters(D,centers,cluster_assignment,n,m,K,euclidean_distance_tmp);

        // 2. Update Centroids
        lloyd_update_centers(D,centers,cluster_assignment,n,m,K,initial_centers);
        num_iter++;
    }

    ///////////////////////////////////
    //// STEP 3 - Set Return Values
    ///////////////////////////////////

    PutRNGstate();
    memcpy(centers_r,centers,K*m*sizeof(double));
    memcpy(cluster_assignment_r,cluster_assignment,n*sizeof(int));
    clock_gettime(CLOCK_MONOTONIC, &cur_time);
    run_time = time_difference_ms(&start_time,&cur_time);
    *time_R = run_time;
    *num_iters_R = num_iter;
}

// REQUIRES iwork of length K
void lloyd_update_centers(double const* restrict D, double* restrict centers, int* restrict cluster_assignment, 
                            int n, int m, int K, int* restrict iwork){
    int tmp1;
    double dtmp1;
    double* restrict tmp_ptr1;
    double const* restrict tmp_ptr2;
    int* restrict cluster_sizes;
    cluster_sizes = iwork;

    // zero out current cluster centers, sizes
    for(int i=0; i < m*K; i++) centers[i] = 0;
    for(int i=0; i < K; i++) cluster_sizes[i] = 0;

    // create new centers
    for(int i=0; i < n; i++){
        tmp1 = cluster_assignment[i];
        cluster_sizes[tmp1] = cluster_sizes[tmp1] + 1;
        tmp_ptr1 = centers + tmp1*m;
        tmp_ptr2 = D + i*m;
        for(int j = 0; j < m; j++){
            dtmp1 = *tmp_ptr1 + *tmp_ptr2;
            *tmp_ptr1 = dtmp1;
            tmp_ptr1++;
            tmp_ptr2++;
        }
    }

    for(int k=0; k < K; k++){
        for(int i=0; i < m; i++){
            tmp1 = k*m + i;
            dtmp1 = centers[tmp1];
            centers[tmp1] = dtmp1 / cluster_sizes[k];
        }
    }
}

int lloyd_update_clusters(double const* restrict D, double* restrict centers, int* restrict cluster_assignment,
                            int n, int m, int K, double* restrict euclidean_distance_tmp){
    double dtmp1;
    int unchanged = 1;

    #pragma omp parallel for shared(unchanged)
    for(int i=0; i < n; i++) {
        int k=0;
        double min_dist = euclidean_distance(D + i*m,centers + k*m,m,euclidean_distance_tmp);
        int min_dist_idx = k;
        for(k = 1; k < K; k++){
            dtmp1 = euclidean_distance(D + i*m,centers + k*m,m,euclidean_distance_tmp);
            if(dtmp1 < min_dist){
                min_dist = dtmp1;
                min_dist_idx = k;
            }
        }
        if(cluster_assignment[i] != min_dist_idx){
            cluster_assignment[i] = min_dist_idx;
            unchanged = 0;
        }
    }
    return unchanged;
}

// MUST BE CALLED FROM WITHIN GetRNGState, PutRNGState PAIR
int sample_discrete_distribution(double* restrict prob_dist,int n){
    double q = unif_rand();
    int cur_idx = 0;
    double total = 0.0;
    for(cur_idx = 0; cur_idx < n; cur_idx++){
        total = total + prob_dist[cur_idx];
        if(total > q) break;
    }
    return cur_idx;
}

void update_min_distance(double const* restrict D, double* restrict min_center_distance,
                            int new_center_idx, int n, int m,double* restrict euclidean_distance_tmp) {
    double dtmp1;

    #pragma omp parallel for
    for(int i=0; i < n; i++){
        //only update if non-zero
        if(min_center_distance[i] > 0.0){
            dtmp1 = euclidean_distance(D+i*m,D+new_center_idx*m,m,euclidean_distance_tmp);
            min_center_distance[i] = min_center_distance[i] > dtmp1 ? dtmp1 : min_center_distance[i];
        }
    }
}

void min_distance_to_probability(double* restrict min_distances, double* restrict prob_dist, int n) {
    double tmp1;
    double dist_sum = 0.0;

    for(int i = 0; i < n; i++){
        dist_sum = dist_sum + min_distances[i];
    }

    #pragma omp simd
    for(int i=0; i < n; i++){
        tmp1 = min_distances[i] / dist_sum;
        prob_dist[i] = tmp1;
    }
}

