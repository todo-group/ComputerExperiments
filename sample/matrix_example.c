#include "cmatrix.h"
#include <stdio.h>

int main(void)
{
    const int m = 2;
    const int n = 3;

    /* 倍精度行列を動的に確保 */
    double** mat = alloc_dmatrix(m, n);
    mat_elem(mat, 0, 0) = 2;
    mat_elem(mat, 0, 1) = 2.5;
    mat_elem(mat, 0, 2) = 1.5;
    mat_elem(mat, 1, 0) = 1;
    mat_elem(mat, 1, 1) = 2;
    mat_elem(mat, 1, 2) = 0.5;

    /* 長さnの倍精度ベクトルを動的に確保 */
    double* v1 = alloc_dvector(n);
    v1[0] = 0.5;
    v1[1] = 1.5;
    v1[2] = 1;

    /* 長さmの倍精度ベクトルを動的に確保 */
    double* v2 = alloc_dvector(m);

    /* v2 = mat * v1 を計算 */
    for (int i = 0; i < m; ++i)
    {
        v2[i] = 0;
        for (int j = 0; j < n; ++j)
        {
            v2[i] += mat_elem(mat, i, j) * v1[j];
        }
    }

    /* 結果を出力 */
    printf("mat = ");
    fprint_dmatrix(stdout, m, n, mat);
    printf("v1 = ");
    fprint_dvector(stdout, n, v1);
    printf("v2 = ");
    fprint_dvector(stdout, m, v2);

    /* 確保したメモリを開放 */
    free_dmatrix(mat);
    free_dvector(v1);
    free_dvector(v2);
    return 0;
}
