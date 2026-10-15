#include "cmatrix.h"
#include <stdio.h>

/* 正方格子の構造を表す配列の生成
   neighbor: (L*L) x 4 の配列
     mat_elem(neighbor, i, k) [k=0,1,2,3] はサイトiの4つの
     最近接サイトを表す
   bond: 2*(L*L) x 2 の配列
     mat_elem(bond, j, 0) と mat_elem(bond, j, 1) はボンドjの
     両端のサイトを表す */

void init_square_lattice(int L, int** neighbor, int** bond)
{
    int b = 0;
    for (int y = 0; y < L; ++y)
    {
        for (int x = 0; x < L; ++x)
        {
            int i, xn, yn, j;
            i = L * y + x;
            /* right */
            xn = (x + 1) % L;
            yn = y;
            j = L * yn + xn;
            mat_elem(neighbor, i, 0) = j;
            mat_elem(bond, b, 0) = i;
            mat_elem(bond, b, 1) = j;
            ++b;
            /* left */
            xn = (L + x - 1) % L;
            yn = y;
            j = L * yn + xn;
            mat_elem(neighbor, i, 1) = j;
            /* up */
            xn = x;
            yn = (y + 1) % L;
            j = L * yn + xn;
            mat_elem(neighbor, i, 2) = j;
            mat_elem(bond, b, 0) = i;
            mat_elem(bond, b, 1) = j;
            ++b;
            /* down */
            xn = x;
            yn = (L + y - 1) % L;
            j = L * yn + xn;
            mat_elem(neighbor, i, 3) = j;
        }
    }
}

/* リストの出力 */

int main(void)
{
    const int L = 4;
    int** neighbor = alloc_imatrix(L * L, 4);
    int** bond = alloc_imatrix(2 * L * L, 2);
    init_square_lattice(L, neighbor, bond);
    printf("L = %d:\n", L);
    printf("neighbor list:\n");
    for (int i = 0; i < L * L; ++i)
    {
        const int x = i % L;
        const int y = i / L;
        printf("  %d (%d,%d): %d %d %d %d\n", i, x, y,
            mat_elem(neighbor, i, 0), mat_elem(neighbor, i, 1),
            mat_elem(neighbor, i, 2), mat_elem(neighbor, i, 3));
    }
    printf("bond list:\n");
    for (int b = 0; b < 2 * L * L; ++b)
    {
        printf("  %d: %d --- %d \n", b,
            mat_elem(bond, b, 0), mat_elem(bond, b, 1));
    }
    return 0;
}
