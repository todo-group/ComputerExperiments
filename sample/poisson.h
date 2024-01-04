#include "cmatrix.h"
#include <math.h>

/* 二次元ポアソン方程式の行列表示: Ax=b

  * 各辺をn等分する (メッシュ幅 h=1/n)
  * (n+1) x (n+1) の二次元メッシュ
  * xとyはそれぞれ0〜nの値を取る
  * x=0,n あるいは y=0,n (合計4n個)は境界条件で決まる固定の値を持つ
    (境界条件は右辺のベクトルbに含める)
  * 残り(n-1)x(n-1)個を一列にならべてベクトルを作る
  * メッシュ座標(x,y)とベクトルのインデックスiの関係
    i = (n-1) * (y-1) + (x-1)
*/

/* 二次元ポワソン方程式の行列次元 */
int matrix_dimension(int n) {
  return (n-1) * (n-1);
}
  
/* 境界条件を与える関数 */
double boundary(int n, int x, int y) {
  if (x == 0) {
    return sin(M_PI * y / n);
  } else if (x == n) {
    return 0;
  } else if (y == 0) {
    return 0;
  } else if (y == n) {
    return 0;
  }
  return 0;
}

/* ポアソン方程式の右辺(電荷密度など)を与える関数 */
double f(int n, int x, int y) {
  /* 例: f=0 (ラプラス方程式) */
  return 0;
}

/* 二次元ポワソン方程式の行列を生成する */
/* 行列: [(n-1)*(n-1)] x [(n-1)*(n-1)] 行列 */
void generate_dense(int n, double h, double **a) {
  int dim = matrix_dimension(n);
  int i, j, x, y, xn, yn;
  double h2inv = 1.0 / (h*h);

  /* ゼロで初期化 */
  for (j = 0; j < dim; ++j) {
    for (i = 0; i < dim; ++i) {
      mat_elem(a, i, j) = 0.0;
    }
  }

  /* 非ゼロの要素を設定 */
  for (x = 1; x < n; ++x) {
    for (y = 1; y < n; ++y) {
      i = (n-1) * (y-1) + (x-1);
      /* diagonal element */
      mat_elem(a, i, i) = -4 * h2inv;
      /* up direction */
      xn = x;
      yn = y + 1;
      if (yn != n) mat_elem(a, i, (n-1) * (yn-1) + (xn-1)) = h2inv;
      /* down direction */
      xn = x;
      yn = y - 1;
      if (yn != 0) mat_elem(a, i, (n-1) * (yn-1) + (xn-1)) = h2inv;
      /* right direction */
      xn = x + 1;
      yn = y;
      if (xn != n) mat_elem(a, i, (n-1) * (yn-1) + (xn-1)) = h2inv;
      /* left direction */
      xn = x - 1;
      yn = y;
      if (xn != 0) mat_elem(a, i, (n-1) * (yn-1) + (xn-1)) = h2inv;
    }
  }
}

/* 二次元ポワソン方程式の行列を生成する */
/* 非ゼロの要素のみを格納
     n_elem[n-1]: 各行の非ゼロ要素の数(最大5)
     col[n-1, 5]: 非ゼロ要素の列
     val[n-1, 5]: 非ゼロ要素の値 */
void generate_sparse(int n, double h, int *n_elem, int **col, double **val) {
  int i, x, y, xn, yn;
  double h2inv = 1.0 / (h*h);
  for (x = 1; x < n; ++x) {
    for (y = 1; y < n; ++y) {
      i = (n-1) * (y-1) + (x-1);
      n_elem[i] = 0;
      /* diagonal element */
      mat_elem(col, i, n_elem[i]) = i;
      mat_elem(val, i, n_elem[i]) = -4 * h2inv;
      n_elem[i] += 1;
      /* up direction */
      xn = x;
      yn = y + 1;
      if (yn != n) {
        mat_elem(col, i, n_elem[i]) = (n-1) * (yn-1) + (xn-1);
        mat_elem(val, i, n_elem[i]) = h2inv;
        n_elem[i] += 1;
      }
      /* down direction */
      xn = x;
      yn = y - 1;
      if (yn != 0)  {
        mat_elem(col, i, n_elem[i]) = (n-1) * (yn-1) + (xn-1);
        mat_elem(val, i, n_elem[i]) = h2inv;
        n_elem[i] += 1;
      }
      /* right direction */
      xn = x + 1;
      yn = y;
      if (xn != n)  {
        mat_elem(col, i, n_elem[i]) = (n-1) * (yn-1) + (xn-1);
        mat_elem(val, i, n_elem[i]) = h2inv;
        n_elem[i] += 1;
      }
      /* left direction */
      xn = x - 1;
      yn = y;
      if (xn != 0)  {
        mat_elem(col, i, n_elem[i]) = (n-1) * (yn-1) + (xn-1);
        mat_elem(val, i, n_elem[i]) = h2inv;
        n_elem[i] += 1;
      }
    }
  }
}

/* 二次元ポワソン方程式の右辺のベクトルを生成する */
/* ベクトル: b[(n-1)*(n-1)] */
void generate_rhs(int n, double h, double *b) {
  int dim = matrix_dimension(n);
  int i, x, y, xn, yn;
  double h2inv = 1.0 / (h*h);

  /* ゼロで初期化 */
  for (i = 0; i < dim; ++i) b[i] = 0.0;

  /* 非ゼロの要素を設定 */
  for (x = 1; x < n; ++x) {
    for (y = 1; y < n; ++y) {
      i = (n-1) * (y-1) + (x-1);
      /* f term */
      b[i] = f(n, x, y);
      /* up direction */
      xn = x;
      yn = y + 1;
      if (yn == n) b[i] -= h2inv * boundary(n, xn, yn);
      /* down direction */
      xn = x;
      yn = y - 1;
      if (yn == 0) b[i] -= h2inv * boundary(n, xn, yn);
      /* right direction */
      xn = x + 1;
      yn = y;
      if (xn == n) b[i] -= h2inv * boundary(n, xn, yn);
      /* left direction */
      xn = x - 1;
      yn = y;
      if (xn == 0) b[i] -= h2inv * boundary(n, xn, yn);
    }
  }
}

void product_sparse(int n, int *n_elem, int **col, double **val,
                    double *v_in, double *v_out) {
  int dim = matrix_dimension(n);
  int i, j;
  for (i = 0; i < dim; ++i) {
    v_out[i] = 0;
    for (j = 0; j < n_elem[i]; ++j) {
      v_out[i] += mat_elem(val, i, j) * v_in[mat_elem(col, i, j)];
    }
  }
}
