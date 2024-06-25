#include "cmatrix.h"
#include "dgesvd.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int imin(int x, int y) { return (x < y) ? x : y; }
int imax(int x, int y) { return (x > y) ? x : y; }

int main(int argc, char** argv)
{

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
        exit(1);
    }

    /* read matrix A from a file */
    FILE* fp = fopen(argv[1], "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: file can not open\n");
        exit(1);
    }
    int m, n;
    double** a;
    read_dmatrix(fp, &m, &n, &a);
    printf("Matrix A:\n");
    fprint_dmatrix(stdout, m, n, a);

    /* allocate matrices and vectors */
    const int r = imin(m, n);
    double** u = alloc_dmatrix(m, r);
    double** vt = alloc_dmatrix(r, n);
    double* s = alloc_dvector(r);
    const int lwork = imax(3 * r + imax(m, n), 5 * r);
    double* work = alloc_dvector(lwork);

    /* perform SVD */
    int info;
    const char jobu = 'S';
    const char jobvt = 'S';
    dgesvd_(&jobu, &jobvt, &m, &n, mat_ptr(a), &m, vec_ptr(s), mat_ptr(u), &m,
        mat_ptr(vt), &r, vec_ptr(work), &lwork, &info);
    if (info != 0)
    {
        fprintf(stderr, "Error: LAPACK::dgesvd failed\n");
        exit(1);
    }
    printf("Result of SVD U:\n");
    fprint_dmatrix(stdout, m, r, u);
    printf("Result of SVD S:\n");
    fprint_dvector(stdout, r, s);
    printf("Result of SVD Vt:\n");
    fprint_dmatrix(stdout, r, n, vt);

    // check the result of SVD
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            mat_elem(a, i, j) = 0.0;
            for (int k = 0; k < r; ++k)
            {
                mat_elem(a, i, j) += mat_elem(u, i, k) * s[k] * mat_elem(vt, k, j);
            }
        }
    }
    printf("Reconstruction of the original matrix A:\n");
    fprint_dmatrix(stdout, m, n, a);

    // approximate A by rank (r-1) matrix
    s[r - 1] = 0.0; // set the last singular value to zero
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            mat_elem(a, i, j) = 0.0;
            for (int k = 0; k < r; ++k)
            {
                mat_elem(a, i, j) += mat_elem(u, i, k) * s[k] * mat_elem(vt, k, j);
            }
        }
    }
    printf("Rank (r-1) approximation of A:\n");
    fprint_dmatrix(stdout, m, n, a);

    free_dmatrix(a);
    free_dmatrix(u);
    free_dmatrix(vt);
    free_dvector(s);
    free_dvector(work);
    return 0;
}
