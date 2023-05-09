#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/* 長さnの配列vの和を求める */
double sum(int n, const double* v)
{
    double s = 0;
    for (int i = 0; i < n; ++i)
    {
        s += v[i];
    }
    return s;
}

int main(void)
{
    const int n = 10;

    /* 長さ10の動的配列 */
    double* v = (double*)malloc((size_t)(n * sizeof(double)));
    assert(v != NULL);

    /* v[0]〜v[9]に1〜10を代入 (v[10]は存在しないことに注意) */
    for (int i = 0; i < n; ++i)
    {
        v[i] = i + 1;
    }

    /* vを表示 */
    printf("v: ");
    for (int i = 0; i < n; ++i)
    {
        printf("%f ", v[i]);
    }
    printf("\n");

    double s;

    /* vの要素の和 (先頭のアドレスを渡す) */
    s = sum(n, v);
    printf("sum of v[0]...v[9] = %f\n", s);

    /* &v[0]と書いても同じ */
    s = sum(n, &v[0]);
    printf("sum of v[0]...v[9] = %f\n", s);

    /* v[1]〜v[9]の和 (v[1]のアドレスを渡す) */
    s = sum(n - 1, v + 1);
    printf("sum of v[1]...v[9] = %f\n", s);

    /* あるいは */
    s = sum(n - 1, &v[1]);
    printf("sum of v[1]...v[9] = %f\n", s);

    /* v[j]〜v[9]の和はどう書ける? */

    /* 動的に確保したメモリ領域の解放 */
    free(v);
    return 0;
}
