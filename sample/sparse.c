/* Hamiltonian of one-dimensional transverse-field Ising model */

/* store only diagonal and non-zero offdiagonal elements */

#include "cmatrix.h"
#include <stdio.h>
#include <stdlib.h>

/* calculate spin state at j-th site from bit pattern c */
int spin(int c, int j)
{
    return 1 - 2 * ((c >> j) & 1);
}

/* flip j-th spin */
int flip(int c, int j)
{
    return c ^ (1 << j);
}

/* fill diag[dim], od_num[dim], od_index[dim][nmax], od_value[dim][nmax] */
void set_matrix(int n, double J, double h, double g, double* diag, int* od_num,
    int** od_index, double** od_value)
{
    const int dim = (1 << n);
    for (int c = 0; c < dim; ++c)
    {
        diag[c] = 0;
        od_num[c] = 0;
    }
    for (int c = 0; c < dim; ++c)
    {
        for (int i = 0; i < n; ++i)
        {
            const int j = (i + 1) % n; /* periodic boundary condition */
            const int s0 = spin(c, i);
            const int s1 = spin(c, j);
            diag[c] += -J * s0 * s1 - h * s0;
        }
        for (int i = 0; i < n; ++i)
        {
            mat_elem(od_index, c, od_num[c]) = flip(c, i);
            mat_elem(od_value, c, od_num[c]) = -g;
            od_num[c] += 1;
        }
    }
}

void multiply(int dim, const double* diag, const int* od_num, int** od_index,
    double** od_value, const double* v, double* w)
{
    for (int c = 0; c < dim; ++c)
    {
        w[c] = diag[c] * v[c];
        for (int i = 0; i < od_num[c]; ++i)
        {
            w[c] += mat_elem(od_value, c, i) * v[mat_elem(od_index, c, i)];
        }
    }
}

int main(int argc, char** argv)
{
    /* default value */
    int n = 2;          /* number of spins */
    const double J = 1; /* coupling constant */
    double h = 0;       /* longitudial field */
    double g = 0;       /* transverse field */

    if (argc > 1)
    {
        n = atoi(argv[1]);
    }
    if (argc > 2)
    {
        h = atof(argv[2]);
    }
    if (argc > 3)
    {
        g = atof(argv[3]);
    }
    const int dim = 1 << n; /* dimension of Hamiltonian */
    const int nmax = n;     /* maxmum number of non-zero offdiagonal elements for each row */
    printf("n = %d; dim = %d; J = %f; h = %f; g = %f\n", n, dim, J, h, g);

    double* diag = alloc_dvector(dim);
    int* od_num = alloc_ivector(dim);
    int** od_index = alloc_imatrix(dim, nmax);
    double** od_value = alloc_dmatrix(dim, nmax);
    set_matrix(n, J, h, g, diag, od_num, od_index, od_value);

    printf("calculate product of Hamiltonian and identity matrix:\n");
    double* v = alloc_dvector(dim);
    double* w = alloc_dvector(dim);
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            v[j] = 0;
        }
        v[i] = 1;
        multiply(dim, diag, od_num, od_index, od_value, v, w);
        fprint_dvector(stdout, dim, w);
    }
    free_dvector(v);
    free_dvector(w);

    free_dvector(diag);
    free_ivector(od_num);
    free_imatrix(od_index);
    free_dmatrix(od_value);
    return 0;
}
