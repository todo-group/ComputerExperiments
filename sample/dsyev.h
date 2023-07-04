#ifndef DSYEV_H
#define DSYEV_H

/* http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html */
void dsyev_(const char* JOBZ, const char* UPLO, const int* N, double* A, const int* LDA, double* W,
    double* WORK, const int* LWORK, int* INFO);

#endif // DSYEV_H