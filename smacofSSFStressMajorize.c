#include "smacofSSFStress.h"

void smacofSSFStressMajorize(int* nobj, int* ndim, int* ndat, double* snew,
                             int* iind, int* jind, double* baux, double* vinv,
                             double* dhat, double* dold, double* dnew,
                             double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k], js = jind[k];
        double elem = baux[k];
        for (int s = 0; s < Ndim; s++) {
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
            is += Nobj;
            js += Nobj;
        }
    }
    int k = 0;
    for (int j = 0; j < Nobj - 1; j++) {
        for (int i = j + 1; i < Nobj; i++) {
            double elem = vinv[k];
            int is = i, js = j;
            for (int s = 0; s < Ndim; s++) {
                double add = elem * (xtmp[is] - xtmp[js]);
                xnew[is] += add;
                xnew[js] -= add;
                is += Nobj;
                js += Nobj;
            }
            k++;
        }
    }
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k], js = jind[k];
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(xnew[is] - xnew[js]);
            is += Nobj;
            js += Nobj;
        }
        dnew[k] = sqrt(sum);
    }
    xfree(xtmp);
    return;
}