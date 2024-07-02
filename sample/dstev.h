#ifndef DSTEV_H
#define DSTEV_H

// LAPACK: DSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices
void dstev_(const char*	JOBZ, const int* N, double* D, double*	E, double* Z, const int* LDZ,
    double* WORK, int* INFO);

#endif // DSTEV_H