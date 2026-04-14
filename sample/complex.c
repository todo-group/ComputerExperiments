#include <complex.h>
#include <math.h>
#include <stdio.h>

int main(void)
{
    const double complex x = 0 + 1 * I; /* 虚数単位 */
    const double complex y = cexp(x * M_PI);
    printf("i = (%lf,%lf)\n", creal(x), cimag(x));
    printf("e^{i*pi} = (%lf,%lf)\n", creal(y), cimag(y));
    return 0;
}
