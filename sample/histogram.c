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
    double* hist = alloc_dvector(bins);
    for (int j = 0; j < bins; ++j)
    {
        hist[j] = 0.0;
    }

    /* generate histogram */
    int n_outside = 0;
    for (int i = 0; i < samples; ++i)
    {
        const double x = ran();
        const int j = (x - xmin) / dx;
        if (x >= xmin && x < xmax && j >= 0 && j < bins)
        {
            hist[j] += 1.0;
        }
        else
        {
            ++n_outside;
        }
    }

    for (int j = 0; j < bins; ++j)
    {
        const double x = xmin + (j + 0.5) * dx;
        const double p = hist[j] / samples;
        const double density = p / dx;
        const double stderr_density = sqrt(p * (1.0 - p) / (samples - 1)) / dx;
        printf("%d %15.10f %15.10f %15.10f\n", j, x, density, stderr_density);
    }

    fprintf(stderr, "outside: %d / %d (%.6f)\n", n_outside, samples, (double)n_outside / samples);

    free_dvector(hist);
    return 0;
}
