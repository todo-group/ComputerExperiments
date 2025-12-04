#include "poisson.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  int n, dim, i, j;
  int *n_elem;
  int **col;
  double **val;
  double *b;
  double *v_in, *v_out;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s n\n", argv[0]);
    exit(1);
  }
  n = atoi(argv[1]);

  dim = matrix_dimension(n);
  n_elem = alloc_ivector(dim);
  col = alloc_imatrix(dim, 5);
  val = alloc_dmatrix(dim, 5);
  generate_sparse(n, 1.0/n, n_elem, col, val);
  b = alloc_dvector(dim);
  generate_rhs(n, 1.0/n, b);

  /* print out matrix */
  v_in = alloc_dvector(dim);
  v_out = alloc_dvector(dim);
  printf("Matrix A:\n");
  printf("%d %d\n", dim, dim);
  for (i = 0; i < dim; ++i) {
    for (j = 0; j < dim; ++j) v_in[j] = 0;
    v_in[i] = 1;
    product_sparse(n, n_elem, col, val, v_in, v_out);
    for (j = 0; j < dim; ++j) printf("%10.5f ", v_out[j]);
    printf("\n");
  }
  printf("Vector B (transposed):\n");
  fprint_dvector(stdout, dim, b);
  
  free_ivector(n_elem);
  free_imatrix(col);
  free_dmatrix(val);
  free_dvector(b);
  free_dvector(v_in);
  free_dvector(v_out);
}
