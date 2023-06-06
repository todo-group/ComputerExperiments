#ifndef DGETRS_H
#define DGETRS_H

/* http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html */
void dgetrs_(const char* TRANS, const int* N, const int* NRHS, double* A, const int* LDA, const int* IPIV, double* B, const int* LDB, int* INFO);

#endif // DGETRS_H