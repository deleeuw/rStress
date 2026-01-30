#ifndef SMACOF_SS_RSTRESS_H
#define SMACOF_SS_RSTRESS_H

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

void smacofMPInverseV(int* nobj, int* ndat, int* iind, int* jind, double* wght,
                      double* vinv);

double smacofSSRStressLoss(int* ndat, double* edis, double* dhat,
                           double* wght, double *rpow);

void smacofSSRStressEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal, int* weighted,
                           double* sold, double* snew, double* eps, double *rpow,
                           int* iind, int* jind, int* blks, double* wght, 
                           double* edis, double* dhat, double* xold, double* xnew);

void smacofSSRStressMajorize(int* nobj, int* ndim, int* ndat, double* snew, int* iind,
                      int* jind, int* weighted, double* wght, double* vinv, double* edis,
                      double* dhat, double* xold, double* xnew);

void smacofSSRStressMonotone(int* ndat, int* ties, double* snew,
                       int* iind, int* jind, int* blks, double* edis,
                       double* dhat, double* wght);

void smacofSSEngine(const int* nobj, const int* ndim, const int* ndat,
                    const int* nord, const int* safe, int* itel, int* kord,
                    const int* ties, const int* itmax, const int* digits,
                    const int* width, const int* verbose, const int* ordinal,
                    const int* weighted, double* sold, double* snew,
                    const double* eps, int* iind, int* jind, int* iord,
                    int* blks, double* wght, double* edis, double* dhat,
                    double* xold, double* xnew);

void smacofSSMajorize(const int* nobj, const int* ndim, const int* ndat,
                      const int* itel, int* kord, const int* nord, int* iind,
                      int* jind, const int* iord, const int* safe,
                      const int* weighted, double* wght, double* vinv,
                      double* dhat, double* xold, double* xnew);

void smacofSSMonotone(const int* ndat, const int* ties, int* iind, int* jind,
                      int* blks, double* edis, double* dhat, double* wght);

void smacofSSSMatrixPrint(double *mat, int *nobj, int *ndat, int *iind, int *jind, int *width, int *digits);


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

#endif /* SMACOF_SS_RSTRESS_H */