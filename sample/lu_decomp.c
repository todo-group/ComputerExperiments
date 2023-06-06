#include "cmatrix.h"
#include "dgetrf.h"
#include "dgetrs.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
        exit(1);
    }

    /* read matrix A and vector b from a file */
    FILE* fp = fopen(argv[1], "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: file can not open\n");
        exit(1);
    }
    int m, n;
    double** matA;
    read_dmatrix(fp, &m, &n, &matA);
    if (m != n)
    {
        fprintf(stderr, "Error: inconsistent number of equations\n");
        exit(1);
    }
    double* b;
    read_dvector(fp, &n, &b);
    if (m != n)
    {
        fprintf(stderr, "Error: inconsistent number of equations\n");
        exit(1);
    }
    printf("Matrix A:\n");
    fprint_dmatrix(stdout, n, n, matA);
    printf("Vector b (transposed):\n");
    fprint_dvector(stdout, n, b);

    /* perform LU decomposition */
    int* ipiv = alloc_ivector(n);
    int info;
    dgetrf_(&n, &n, mat_ptr(matA), &n, vec_ptr(ipiv), &info);
    if (info != 0)
    {
        fprintf(stderr, "Error: LAPACK::dgetrf failed\n");
        exit(1);
    }
    printf("Result of LU decomposition:\n");
    fprint_dmatrix(stdout, n, n, matA);
    printf("Pivot for LU decomposition:\n");
    fprint_ivector(stdout, n, ipiv);

    /* calculate L, U, P */
    double** matL = alloc_dmatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            mat_elem(matL, i, j) = mat_elem(matA, i, j);
        }
        mat_elem(matL, i, i) = 1;
    }
    printf("L:\n");
    fprint_dmatrix(stdout, n, n, matL);

    double** matU = alloc_dmatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            mat_elem(matU, i, j) = mat_elem(matA, i, j);
        }
    }
    printf("U:\n");
    fprint_dmatrix(stdout, n, n, matU);

    double** matP = alloc_dmatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        mat_elem(matP, i, i) = 1;
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const double tp = mat_elem(matP, ipiv[i] - 1, j);
            mat_elem(matP, ipiv[i] - 1, j) = mat_elem(matP, i, j);
            mat_elem(matP, i, j) = tp;
        }
    }
    printf("P:\n");
    fprint_dmatrix(stdout, n, n, matP);

    double** work = alloc_dmatrix(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                mat_elem(work, i, j) += mat_elem(matL, i, k) * mat_elem(matU, k, j);
            }
        }
    }
    printf("L * U:\n");
    fprint_dmatrix(stdout, n, n, work);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            mat_elem(work, i, j) = 0;
            for (int k = 0; k < n; ++k)
            {
                for (int l = 0; l < n; ++l)
                {
                    mat_elem(work, i, j) += mat_elem(matP, k, i) * mat_elem(matL, k, l) * mat_elem(matU, l, j);
                }
            }
        }
    }
    printf("P^t * L * U:\n");
    fprint_dmatrix(stdout, n, n, work);

    free_dmatrix(matL);
    free_dmatrix(matU);
    free_dmatrix(matP);
    free_dmatrix(work);

    /* solve equations */
    const char trans = 'N';
    const int nrhs = 1;
    dgetrs_(&trans, &n, &nrhs, mat_ptr(matA), &n, vec_ptr(ipiv), vec_ptr(b), &n, &info);
    if (info != 0)
    {
        fprintf(stderr, "Error: LAPACK::dgetrs failed\n");
        exit(1);
    }
    printf("Solution X (transposed):\n");
    fprint_dvector(stdout, n, b);

    free_dmatrix(matA);
    free_dvector(b);
    free_ivector(ipiv);
    return 0;
}
