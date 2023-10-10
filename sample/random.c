#include "mersenne_twister.h"
#include <stdio.h>

int main(void)
{
    const int seed = 12345; /* 乱数の種。種が同じ場合には同じ乱数列が生成される */
    const int n = 100;

    /* set seed of Mersenne-Twister pseudo random number generator */
    init_genrand(seed);

    /* generate (0,1) random numbers */
    for (int i = 0; i < n; ++i)
    {
        printf("%d %15.10f\n", i, genrand_real3());
    }
    return 0;
}
