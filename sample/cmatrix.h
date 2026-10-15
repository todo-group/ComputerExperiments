/* Utility for vector and matrix allocation/deallocation */

/* List of functions

  alloc_dvector : allocate vector of double
  alloc_fvector : allocate vector of float
  alloc_zvector : allocate vector of double complex
  alloc_cvector : allocate vector of float complex
  alloc_ivector : allocate vector of int

  alloc_dmatrix : allocate matrix of double
  alloc_fmatrix : allocate matrix of float
  alloc_zmatrix : allocate matrix of double complex
  alloc_cmatrix : allocate matrix of float complex
  alloc_imatrix : allocate matrix of int

  free_dvector : deallocate vector of double
  free_fvector : deallocate vector of float
  free_zvector : deallocate vector of double complex
  free_cvector : deallocate vector of float complex
  free_ivector : deallocate vector of int

  free_dmatrix : deallocate matrix of double
  free_fmatrix : deallocate matrix of float
  free_zmatrix : deallocate matrix of double complex
  free_cmatrix : deallocate matrix of float complex
  free_imatrix : deallocate matrix of int

  fprint_dvector: print out vector of double
  fprint_fvector: print out vector of float
  fprint_zvector: print out vector of double complex
  fprint_cvector: print out vector of float complex
  fprint_ivector: print out vector of int

  fprint_dmatrix: print out matrix of double
  fprint_fmatrix: print out matrix of float
  fprint_zmatrix: print out matrix of double complex
  fprint_cmatrix: print out matrix of float complex
  fprint_imatrix: print out matrix of int

  read_dvector: read vector of double from file
  read_fvector: read vector of float from file
  read_ivector: read vector of int from file

  read_dmatrix: read matrix of double from file
  read_fmatrix: read matrix of float from file
  read_imatrix: read matrix of int from file
*/

#ifndef CMATRIX_H
#define CMATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

/* useful macros */
#define vec_ptr(vec) &(vec)[0]
#define mat_ptr(mat) &(mat)[0][0]
#define mat_elem(mat, i, j) (mat)[(j)][(i)]

/* allocate vector of double */
static inline double *alloc_dvector(int n) {
  double *vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of float */
static inline float *alloc_fvector(int n) {
  float *vec;
  vec = (float*)calloc(n, sizeof(float));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_fvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of double complex */
static inline double complex *alloc_zvector(int n) {
  double complex *vec;
  vec = (double complex*)calloc(n, sizeof(double complex));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_zvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of float complex */
static inline float complex *alloc_cvector(int n) {
  float complex *vec;
  vec = (float complex*)calloc(n, sizeof(float complex));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_cvector\n");
    exit(1);
  }
  return vec;
}

/* allocate vector of int */
static inline int *alloc_ivector(int n) {
  int *vec;
  vec = (int*)calloc(n, sizeof(int));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_ivector\n");
    exit(1);
  }
  return vec;
}

/* allocate m x n column-major matrix of double */
static inline double **alloc_dmatrix(int m, int n) {
  int i;
  double **mat;
  mat = (double**)malloc((size_t)(n * sizeof(double*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  mat[0] = (double*)calloc(m * n, sizeof(double));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n column-major matrix of float */
static inline float **alloc_fmatrix(int m, int n) {
  int i;
  float **mat;
  mat = (float**)malloc((size_t)(n * sizeof(float*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_fmatrix\n");
    exit(1);
  }
  mat[0] = (float*)calloc(m * n, sizeof(float));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_fmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n column-major matrix of double complex */
static inline double complex **alloc_zmatrix(int m, int n) {
  int i;
  double complex **mat;
  mat = (double complex**)malloc((size_t)(n * sizeof(double complex*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  mat[0] = (double complex*)calloc(m * n, sizeof(double complex));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_zmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n column-major matrix of float complex */
static inline float complex **alloc_cmatrix(int m, int n) {
  int i;
  float complex **mat;
  mat = (float complex**)malloc((size_t)(n * sizeof(float complex*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_fmatrix\n");
    exit(1);
  }
  mat[0] = (float complex*)calloc(m * n, sizeof(float complex));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_cmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* allocate m x n column-major matrix of int */
static inline int **alloc_imatrix(int m, int n) {
  int i;
  int **mat;
  mat = (int**)malloc((size_t)(n * sizeof(int*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_imatrix\n");
    exit(1);
  }
  mat[0] = (int*)calloc(m * n, sizeof(int));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_imatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

/* deallocate vector of double */
static inline void free_dvector(double *vec) {
  free(vec);
}

/* deallocate vector of float */
static inline void free_fvector(float *vec) {
  free(vec);
}

/* deallocate vector of double complex */
static inline void free_zvector(double complex *vec) {
  free(vec);
}

/* deallocate vector of float complex */
static inline void free_cvector(float complex *vec) {
  free(vec);
}

/* deallocate vector of int */
static inline void free_ivector(int *vec) {
  free(vec);
}

/* deallocate matrix of double */
static inline void free_dmatrix(double **mat) {
  free(mat[0]);
  free(mat);
}

/* deallocate float matrix of float */
static inline void free_fmatrix(float **mat) {
  free(mat[0]);
  free(mat);
}

/* deallocate matrix of double complex */
static inline void free_zmatrix(double complex **mat) {
  free(mat[0]);
  free(mat);
}

/* deallocate float matrix of float complex */
static inline void free_cmatrix(float complex **mat) {
  free(mat[0]);
  free(mat);
}

/* deallocate float matrix of int */
static inline void free_imatrix(int **mat) {
  free(mat[0]);
  free(mat);
}

/* print out vector of double */
static inline void fprint_dvector(FILE *fp, int n, double *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, "%10.5f ", vec[i]);
  fprintf(fp, "\n");
}

/* print out vector of float */
static inline void fprint_fvector(FILE *fp, int n, float *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, "%10.5f ", vec[i]);
  fprintf(fp, "\n");
}

/* print out vector of double complex */
static inline void fprint_zvector(FILE *fp, int n, double complex *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, "(%10.5f,%10.5f) ", creal(vec[i]), cimag(vec[i]));
  fprintf(fp, "\n");
}

/* print out vector of float complex */
static inline void fprint_cvector(FILE *fp, int n, float complex *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, "(%10.5f,%10.5f) ", creal(vec[i]), cimag(vec[i]));
  fprintf(fp, "\n");
}

/* print out vector of int */
static inline void fprint_ivector(FILE *fp, int n, int *vec) {
  int i;
  fprintf(fp, "%d\n", n);
  for (i = 0; i < n; ++i) fprintf(fp, "%d ", vec[i]);
  fprintf(fp, "\n");
}

/* print out matrix of double */
static inline void fprint_dmatrix(FILE *fp, int m, int n, double **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, "%10.5f ", mat_elem(mat, i, j));
    fprintf(fp, "\n");
  }
}

/* print out matrix of float */
static inline void fprint_fmatrix(FILE *fp, int m, int n, float **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, "%10.5f ", mat_elem(mat, i, j));
    fprintf(fp, "\n");
  }
}

/* print out matrix of double complex */
static inline void fprint_zmatrix(FILE *fp, int m, int n, double complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, "(%10.5f,%10.5f) ", creal(mat_elem(mat, i, j)), cimag(mat_elem(mat, i, j)));
    fprintf(fp, "\n");
  }
}

/* print out matrix of float complex */
static inline void fprint_cmatrix(FILE *fp, int m, int n, float complex **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, "(%10.5f,%10.5f) ", creal(mat_elem(mat, i, j)), cimag(mat_elem(mat, i, j)));
    fprintf(fp, "\n");
  }
}

/* print out matrix of int */
static inline void fprint_imatrix(FILE *fp, int m, int n, int **mat) {
  int i, j;
  fprintf(fp, "%d %d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) fprintf(fp, "%d ", mat_elem(mat, i, j));
    fprintf(fp, "\n");
  }
}

/* read vector of double from file */
static inline void read_dvector(FILE *fp, int *n, double **vec) {
  int i;
  fscanf(fp, "%d", n);
  *vec = alloc_dvector(*n);
  for (i = 0; i < *n; ++i) fscanf(fp, "%lf", &(*vec)[i]);
}

/* read vector of float from file */
static inline void read_fvector(FILE *fp, int *n, float **vec) {
  int i;
  fscanf(fp, "%d", n);
  *vec = alloc_fvector(*n);
  for (i = 0; i < *n; ++i) fscanf(fp, "%f", &(*vec)[i]);
}

/* read matrix of double from file */
static inline void read_dmatrix(FILE *fp, int *m, int *n, double ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_dmatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%lf", &mat_elem(*mat, i, j));
  }
}

/* read matrix of float from file */
static inline void read_fmatrix(FILE *fp, int *m, int *n, float ***mat) {
  int i, j;
  fscanf(fp, "%d", m);
  fscanf(fp, "%d", n);
  *mat = alloc_fmatrix(*m, *n);
  for (i = 0; i < *m; ++i) {
    for (j = 0; j < *n; ++j) fscanf(fp, "%f", &mat_elem(*mat, i, j));
  }
}

#endif
