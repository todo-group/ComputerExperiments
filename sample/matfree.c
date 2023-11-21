/* Hamiltonian of one-dimensional transverse-field Ising model */

/* calculate elements on the fly */

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

void multiply(int n, double J, double h, double g, const double* v, double* w)
{
    const int dim = (1 << n);
    for (int c = 0; c < dim; ++c)
    {
        double diag = 0;
        for (int i = 0; i < n; ++i)
        {
            const int j = (i + 1) % n; /* periodic boundary condition */
            const int s0 = spin(c, i);
            const int s1 = spin(c, j);
            diag += -J * s0 * s1 - h * s0;
        }
        w[c] = diag * v[c];
        for (int i = 0; i < n; ++i)
        {
            w[c] += -g * v[flip(c, i)];
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
    printf("n = %d; dim = %d; J = %f; h = %f; g = %f\n", n, dim, J, h, g);

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
        multiply(n, J, h, g, v, w);
        fprint_dvector(stdout, dim, w);
    }
    free_dvector(v);
    free_dvector(w);
    return 0;
}
