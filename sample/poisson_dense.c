#include "poisson.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  int n, dim;
  double **a;
  double *b;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s n\n", argv[0]);
    exit(1);
  }
  n = atoi(argv[1]);

  dim = matrix_dimension(n);
  a = alloc_dmatrix(dim, dim);
  generate_dense(n, 1.0/n, a);
  b = alloc_dvector(dim);
  generate_rhs(n, 1.0/n, b);

  printf("Matrix A:\n");
  fprint_dmatrix(stdout, dim, dim, a);
  printf("Vector B (transposed):\n");
  fprint_dvector(stdout, dim, b);
  
  free_dmatrix(a);
  free_dvector(b);
}
