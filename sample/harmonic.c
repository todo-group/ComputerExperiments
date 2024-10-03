#include "mersenne_twister.h"
#include <math.h>
#include <stdio.h>

double v(double x)
{
    return x * x;
}

int main(void)
{
    const int seed = 12345;
    const double beta = 1.0;
    const double delta = 10.0;
    const int n = 1000;

    /* set seed of Mersenne-Twister pseudo random number generator */
    init_genrand(seed);

    /* Monte Carlo steps */
    double x = 0.0;
    double x2 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double trial = x + delta * (2 * genrand_real3() - 1);
        if (genrand_real3() < exp(-beta * (v(trial) - v(x))))
        {
            x = trial;
        }
        x2 += x * x;
        printf("%d %10.5f\n", i, x);
    }

    x2 /= n;
    printf("# average of x^2 = %10.5f\n", x2);
    return 0;
}
