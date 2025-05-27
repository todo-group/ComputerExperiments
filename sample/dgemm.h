#ifndef DGEMM_H
#define DGEMM_H

/* http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html */
void dgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K,
    const double* ALPHA, const double* A, const int* LDA, const double* B, const int* LDB,
    const double* BETA, double* C, const int* LDC);

#endif // DGEMM_H