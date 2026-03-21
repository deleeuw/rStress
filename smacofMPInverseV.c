

#include "smacofSSFStress.h"

/*

int nobj = 3, ndat = 3;
int iind[3] = {2, 1, 2};
int jind[3] = {0, 0, 1};
int width = 10;
int digits = 6;
double wght[3] = {3.0, 2.0, 3.0};
double vinv[3] = {0};

int main(void) {
    (void)smacofMPInverseV(&nobj, &ndat, iind, jind, wght, vinv);
    (void)smacofSSSMatrixPrint(wght, &nobj, &ndat, iind, jind, &width, &digits);
    (void)smacofSSSMatrixPrint(vinv, &nobj, &ndat, iind, jind, &width, &digits);
    return EXIT_SUCCESS;
}

*/

/*
clang -o tester smacofMPInverseV.c smacofSSPrint.c syminv.c cholesky.c
timestamp.c
*/

void smacofMPInverseV(const int* nobj, const int* ndat, const int* iind, const int* jind, const double* wght,
                      double* vinv) {
    int Nobj = *nobj, Ndat = *ndat;
    int Nvec = Nobj * (Nobj + 1) / 2;
    int fault = 0, nullty = 0;
    double add = 1.0 / Nobj;
    double* vvec = xcalloc(Nvec, sizeof(double));
    double* vwrk = xcalloc(Nvec, sizeof(double));
    double* work = xcalloc(Nobj, sizeof(double));
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        int ij = i * (i + 1) / 2 + j;
        vvec[ij] = -wght[k];
    }
    for (int i = 0; i < Nobj; i++) {
        int ii = i * (i + 3) / 2;
        for (int j = 0; j < Nobj; j++) {
            if (j == i) {
                continue;
            }
            int ma = MAX(i, j);
            int mi = MIN(i, j);
            int ij = ma * (ma + 1) / 2 + mi;
            vvec[ii] -= vvec[ij];
        }
    }
    for (int i = 0; i < Nvec; i++) {
        vvec[i] += add;
    }
    (void)syminv(vvec, Nobj, vwrk, work, &nullty, &fault);
    for (int i = 0; i < Nvec; i++) {
        vwrk[i] -= add;
    }
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        int ij = i * (i + 1) / 2 + j;
        vinv[k] = vwrk[ij];
    }
    xfree(vvec);
    xfree(vwrk);
    xfree(work);
}