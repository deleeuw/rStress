#ifndef SMACOF_SS_FSTRESS_H
#define SMACOF_SS_FSTRESS_H

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define EPS 1e-6
#define true 1
#define false 0
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


void primaryApproach(const int *, const int *, double *, double *,
                     double *, int *, int *);
void secondaryApproach(const int *, const int *, double *, double *);
void tertiaryApproach(const int *, const int *, double *, double *);
void tieBlockAverages(const int *, const int *, const int *,
                      const double *, const double *, double *, double *,
                      int *, double *);
void monotone(const int *, double *, double *);
int myComp(const void *, const void *);
void mySort(double *, double *, double *, int *, int *, const int *);

void smacofSSFStressFList(double *x, double *f, double *g, int *what, double* rpow, int *nvec);

void smacofMPInverseV(int* nobj, int* ndat, int* iind, int* jind, double* wght,
                      double* vinv);

void smacofSSFStressEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal,
                           double* sold, double* snew, double* eps, int *what,
                           double* rpow, int* iind, int* jind, int* blks, double* wght, 
                           double* edis, double* dhat, double* xold, double* xnew);

void smacofSSFStressMajorize(int* nobj, int* ndim, int* ndat, double* snew,
                             int* iind, int* jind, double* baux, double* vinv,
                             double* dhat, double* dold, double* dnew,
                             double* xold, double* xnew);

void smacofSSFStressMonotone(int* ndat, int* ties, double* snew,
                       int* iind, int* jind, int* blks, double* edis,
                       double* dhat, double* wght);

void smacofSSSMatrixPrint(double *mat, int *nobj, int *ndat, 
  int *iind, int *jind, int *width, int *digits);

void smacofSSRMatrixPrint(double *mat, int* nrow, int* ncol, 
  int* width, int* digits);

void smacofSSVectorPrint(double* vec, int* n, int* width, int* digits);

void cholesky ( double a[], int n, int nn, double u[], int *nullty, 
  int *ifault );

void syminv ( double a[], int n, double c[], double w[], int *nullty, 
  int *ifault );

void timestamp (void);


static inline void *xmalloc(const size_t size) {
  void *p = malloc(size);
  if (!p && size) {
    fprintf(stderr, "FATAL: malloc(%zu) failed\n", size);
    abort();
  }
  return p;
}

static inline void *xcalloc(const size_t nmemb, const size_t size) {
  if (size && nmemb > SIZE_MAX / size) {
    fprintf(stderr, "FATAL: calloc overflow (%zu,%zu)\n", nmemb, size);
    abort();
  }
  void *p = calloc(nmemb, size);
  if (!p && nmemb && size) {
    fprintf(stderr, "FATAL: calloc(%zu,%zu) failed\n", nmemb, size);
    abort();
  }
  return p;
}

static inline void *xrealloc(void *ptr, const size_t size) {
  void *p = realloc(ptr, size);
  if (!p && size != 0) {
    fprintf(stderr, "FATAL: realloc(%p,%zu) failed\n", ptr, size);
    abort();
  }
  return p;
}

#define xfree(p)                                                               \
  {                                                                            \
    if ((p) != NULL) {                                                         \
      free(p);                                                                 \
      p = NULL;                                                                \
    }                                                                          \
  }

#endif /* SMACOF_SS_FSTRESS_H */