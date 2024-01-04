#include "func_2d.h"
#include <math.h>
#include <stdio.h>

double optimize_1d(double* x, const double* grad)
{
    const double epsilon = 0.1;
    const double w = (3.0 - sqrt(5)) / 2.0; /* golden ratio */
    double a, b, c;
    double xa[2], xb[2], xc[2];

    a = 0;
    for (int i = 0; i < 2; ++i)
    {
        xa[i] = x[i] + a * grad[i];
    }
    c = -epsilon;
    for (int i = 0; i < 2; ++i)
    {
        xc[i] = x[i] + c * grad[i];
    }

    if (f(xc) > f(xa))
    {
        // bisection method
        while (1)
        {
            b = (a + c) / 2;
            for (int i = 0; i < 2; ++i)
            {
                xb[i] = x[i] + b * grad[i];
            }
            if (f(xb) < f(xa))
            {
                break;
            }
            c = b;
            for (int i = 0; i < 2; ++i)
            {
                xc[i] = x[i] + c * grad[i];
            }
        }
    }
    else
    {
        // extend method
        while (1)
        {
            b = c;
            for (int i = 0; i < 2; ++i)
            {
                xb[i] = x[i] + b * grad[i];
            }
            c -= epsilon;
            for (int i = 0; i < 2; ++i)
            {
                xc[i] = x[i] + c * grad[i];
            }
            if (f(xc) > f(xb))
            {
                break;
            }
            a = b;
            for (int i = 0; i < 2; ++i)
            {
                xa[i] = x[i] + a * grad[i];
            }
        }
    }

    // start golden section
    while (fabs(c - a) > 1.0e-6 || (f(xa) - f(xb) > 1.0e-12) || (f(xc) - f(xb) > 1.0e-12))
    {
        if (fabs(c - b) > fabs(b - a))
        {
            const double p = b + w * (c - b);
            double xp[2];
            for (int i = 0; i < 2; ++i)
            {
                xp[i] = x[i] + p * grad[i];
            }
            if (f(xp) < f(xb))
            {
                a = b;
                for (int i = 0; i < 2; ++i)
                {
                    xa[i] = x[i] + a * grad[i];
                }
                b = p;
                for (int i = 0; i < 2; ++i)
                {
                    xb[i] = x[i] + b * grad[i];
                }
            }
            else
            {
                c = p;
                for (int i = 0; i < 2; ++i)
                {
                    xc[i] = x[i] + c * grad[i];
                }
            }
        }
        else
        {
            const double p = b + w * (a - b);
            double xp[2];
            for (int i = 0; i < 2; ++i)
            {
                xp[i] = x[i] + p * grad[i];
            }
            if (f(xp) < f(xb))
            {
                c = b;
                for (int i = 0; i < 2; ++i)
                {
                    xc[i] = x[i] + c * grad[i];
                }
                b = p;
                for (int i = 0; i < 2; ++i)
                {
                    xb[i] = x[i] + b * grad[i];
                }
            }
            else
            {
                a = p;
                for (int i = 0; i < 2; ++i)
                {
                    xa[i] = x[i] + a * grad[i];
                }
            }
        }
    }
    return b;
}

int main(void)
{
    double x[2] = { -15.0, 3.0 };
    const int n = 10000;

    printf("%d %lf %lf %lf\n", 0, x[0], x[1], f(x));
    for (int i = 1; i < n; ++i)
    {
        double grad[2];
        df(x, grad);
        const double alpha = optimize_1d(x, grad);
        x[0] += alpha * grad[0];
        x[1] += alpha * grad[1];
        const double diff = sqrt((x[0] - minx) * (x[0] - minx) + (x[1] - miny) * (x[1] - miny));
        printf("%d %lf %lf %20.12lf %20.12lf\n", i, x[0], x[1], f(x), diff);
        if (grad[0] * grad[0] + grad[1] * grad[1] < 1.0e-10)
        {
            break;
        }
    }
    return 0;
}
