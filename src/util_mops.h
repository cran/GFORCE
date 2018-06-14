#ifndef UTIL_MOPS_H
#define UTIL_MOPS_H

// VECTOR, MATRIX OPS
// Computes A = A + c, where A is a vector of length d, c is scalar
void daps(double* restrict A,int inc_A, double c, int d);
// Computes A = AB, where A is a square matrix, and B is diagonal if side=='R'
// Computes A = BA, where A is a square matrix, and B is diagonal if side=='L'
// B is the diagonal of matrix B
// column major layout means 'R' should be FASTER
void dsmtd(double* restrict A,double* restrict B, int d, const char side);
// Computes entrywise exponential function of A
void dvexp(double* restrict A, int d);
// Computes sum of all entries in vector A, of length d
double dsumv(double* restrict A, int d);
// Computes trace of square matrix A, which is dxd
double dtrace(double* restrict A, int d);
// Compute row or column sums of dXd matrix A, store result in A_rsums
void dcsum(double* restrict A, int d, double* restrict A_rsums);
// Compute X + Y and store in Z
void dxpyez(int d, double* restrict X, double* restrict Y, double* restrict Z);
// Compute y = a*x + b*y
void daxpby(double a, double* restrict X, double b, double* restrict Y, int d);
// Compute x = tr(A * B(G))
double dabgtp(double* restrict A, int* restrict G, int d, int K,int* iwork,double* dwork);


#endif
