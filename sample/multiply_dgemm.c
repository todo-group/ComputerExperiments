#include "cmatrix.h"
#include "dgemm.h"
#include "mersenne_twister.h"
#include <stdio.h>

int main(void)
{
    const int n = 1000;
    const int seed = 12345;

    /* initialize random number generator */
    init_genrand(seed);

    /* generate random matrices */
    double** matA = alloc_dmatrix(n, n);
    double** matB = alloc_dmatrix(n, n);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            mat_elem(matA, i, j) = genrand_real3();
            mat_elem(matB, i, j) = genrand_real3();
        }
    }

    /* calculate C = A * B */
    double** matC = alloc_dmatrix(n, n);
    const char trans = 'N';
    const double alpha = 1.0;
    const double beta = 0.0;
    dgemm_(&trans, &trans, &n, &n, &n, &alpha, mat_ptr(matA), &n, mat_ptr(matB), &n, &beta, mat_ptr(matC), &n);

    printf("n = %d\n", n);
    printf("matC[0][0] = %15.10f\n", mat_elem(matC, 0, 0));
    printf("matC[0][n-1] = %15.10f\n", mat_elem(matC, 0, n - 1));
    printf("matC[n-1][n-1] = %15.10f\n", mat_elem(matC, n - 1, n - 1));
    return 0;
}
