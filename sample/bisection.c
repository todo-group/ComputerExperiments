#include <math.h>
#include <stdio.h>

double f(double x)
{
    return sin(x);
}

int main(void)
{
    const double delta = 1.0e-1;
    const double epsilon = 1.0e-8;
    double a = 0.3 * M_PI;
    double b = a + delta;
    printf("Find initial enclosure\n");
    while (f(a) * f(b) > 0)
    {
        printf("f(%le) = %le, f(%le) = %le\n", a, f(a), b, f(b));
        b += delta;
    }
    printf("Start bisection search\n");
    while (fabs(b - a) > epsilon)
    {
        printf("f(%le) = %le, f(%le) = %le\n", a, f(a), b, f(b));
        const double x = (a + b) / 2;
        if (fabs(f(x)) < epsilon)
        {
            printf("Solution = %le\n", x);
            break;
        }
        if (f(a) * f(x) < 0)
        {
            b = x;
        }
        else
        {
            a = x;
        }
    }
    return 0;
}
