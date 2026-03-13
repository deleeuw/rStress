#include "smacofSSFStress.h"

void smacofSSSMatrixPrint(double* mat, int* nobj, int* ndat, int* iind,
                          int* jind, int* width, int* digits) {
    int Nobj = *nobj, Ndat = *ndat;
    double* out = xcalloc(Nobj * Nobj, sizeof(double));
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        out[i + (Nobj - 1) * j] = out[j + (Nobj - 1) * i] = mat[k];
    }
    for (int i = 0; i < Nobj; i++) {
        for (int j = 0; j < Nobj; j++) {
            printf(" %*.*f ", *width, *digits, out[i + (Nobj - 1) * j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofSSVectorPrint(double* vec, int* n, int* width, int* digits) {
    for (int i = 0; i < *n; i++) {
        printf(" %*.*f ", *width, *digits, vec[i]);
    }
    printf("\n\n");
    return;
}