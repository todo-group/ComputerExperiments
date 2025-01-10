#ifndef DGETRF_H
#define DGETRF_H

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
void dgetrf_(const int* M,const int* N, double* A,const int* LDA, int* IPIV, int* INFO);

#endif // DGETRF_H