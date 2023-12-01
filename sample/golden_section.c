#include "func_1d.h"
#include <stdio.h>

int main(void)
{
    double x = 2.0;
    double dx = 0.01;
    double a, b, c;

    const double w = (3.0 - sqrt(5)) / 2.0; /* golden ratio */

    /* 初期囲い込み開始 */
    if (f(x - dx) > f(x) && f(x + dx) > f(x))
    {
        a = x - dx;
        b = x;
        c = x + dx;
    }
    else if (f(x - dx) < f(x))
    {
        /* 左へ進む */
        c = x;
        b = x - dx;
        dx = 2 * dx;
        while (f(b - dx) < f(b))
        {
            c = b;
            b = b - dx;
            dx = 2 * dx;
        }
        a = b - dx;
    }
    else
    {
        /* 右へ進む */
        a = x;
        b = x + dx;
        dx = 2 * dx;
        while (f(b + dx) < f(b))
        {
            a = b;
            b = b + dx;
            dx = 2 * dx;
        }
        c = b + dx;
    }
    /* 初期囲い込み完了 */
    printf("a = %15.10lf, f(a) = %15.10lf\n", a, f(a));
    printf("b = %15.10lf, f(b) = %15.10lf\n", b, f(b));
    printf("c = %15.10lf, f(c) = %15.10lf\n", c, f(c));

    /* 黄金分割法 */
    while ((c - a) > 1.0e-6 || (f(a) - f(b) > 1.0e-12) || (f(c) - f(b) > 1.0e-12))
    {
        if ((c - b) > (b - a))
        {
            x = b + w * (c - b);
            if (f(x) < f(b))
            {
                a = b;
                b = x;
            }
            else
            {
                c = x;
            }
        }
        else
        {
            x = b + w * (a - b);
            if (f(x) < f(b))
            {
                c = b;
                b = x;
            }
            else
            {
                a = x;
            }
        }
        printf("a = %15.10lf, b = %15.10lf, c = %15.10lf, f(b) = %15.10lf\n", a, b, c, f(b));
    }
    return 0;
}
