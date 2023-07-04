#ifndef DGESVD_H
#define DGESVD_H

/* http://www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f.html */
void dgesvd_(const char* JOBU, const char* JOBVT, const int* M, const int* N, double* A,
    const int* LDA, double* S, double* U, const int* LDU, double* VT, const int* LDVT,
    double* WORK, const int* LWORK, int* INFO);

#endif // DGESVD_H