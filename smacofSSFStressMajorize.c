
#include "smacofSSFStress.h"

void smacofSSFStressMajorize(int* nobj, int* ndim, int* ndat, double* snew,
                             int* iind, int* jind, double* baux, double* vinv,
                             double* dhat, double* dold, double* dnew,
                             double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = xnew[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k], j = jind[k];
        double elem = baux[k];
        for (int s = 0; s < Ndim; s++) {
            int is = i + Nobj * s, js = j + Nobj * s;
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
        }
    }
    int k = 0;
    for (int i = 0; i < Nobj; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                k++;
                continue;
            }
            double elem = -vinv[k];
            for (int s = 0; s < Ndim; s++) {
                int is = i + Nobj * s, js = j + Nobj * s;
                double add = elem * (xtmp[is] - xtmp[js]);
                xnew[is] += add;
                xnew[js] -= add;
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