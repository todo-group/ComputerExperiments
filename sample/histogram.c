#include "cmatrix.h"
#include "mersenne_twister.h"
#include <math.h>
#include <stdio.h>

/* uniform random number on (0,1) */
double ran(void)
{
    double x = genrand_real3();
    return x;
}

int main(void)
{
    const int seed = 12345;
    const int samples = 1024;
    const int bins = 16;
    const double xmin = 0.0;
    const double xmax = 1.0;

    /* set seed of Mersenne-Twister pseudo random number generator */
    init_genrand(seed);

    /* prepare vector for histogram */
    const double dx = (xmax - xmin) / bins;
    const double dh = 1.0 / (samples * dx);
    double* hist = alloc_dvector(bins);
    for (int j = 0; j < bins; ++j)
    {
        hist[j] = 0.0;
    }

    /* generate histogram */
    for (int i = 0; i < samples; ++i)
    {
        const double x = ran();
        const int j = (x - xmin) / dx;
        hist[j] += 1.0;
    }

    for (int j = 0; j < bins; ++j)
    {
        const double x = xmin + (j + 0.5) * dx;
        const double ave = hist[j] * dh;
        const double err = sqrt(hist[j]) * dh;
        printf("%d %15.10f %15.10f %15.10f\n", j, x, ave, err);
    }

    free_dvector(hist);
    return 0;
}
