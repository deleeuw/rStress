#include "smacofSSFStress.h"

void smacofSSRStressEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal, int* weighted,
                           double* sold, double* snew, double* eps, int* what,
                           double* rpow, int* iind, int* jind, int* blks, double* wght,
                           double* edis, double* dhat, double* xold,
                           double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double Eps = *eps;
    double* waux = xmalloc(Ndat * sizeof(double));
    double* daux = xmalloc(Ndat * sizeof(double));
    double* fval = xmalloc(Ndat * sizeof(double));
    double* gval = xmalloc(Ndat * sizeof(double));
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    while (true) {
        double told = 0.0;
        (void)smacofSSFStressFList(edis, fval, gval, what, rpow, ndat);
        for (int k = 0; k < Ndat; k++) {
            waux[k] = wght[k] * SQUARE(gval[k]);
            daux[k] = (dhat[k] - fval[k]) / gval[k] + edis[k];
            told += waux[k] * SQUARE(daux[k] - edis[k]);
        }
        (void)smacofMPInverseV(nobj, ndat, iind, jind, waux, vinv);
        (void)smacofSSFStressMajorize(nobj, ndim, ndat, snew, iind, jind,
                                      weighted, waux, vinv, edis, daux, xold,
                                      xnew);
        double tnew = 0.0, smid = 0.0;
        for (int k = 0; k < Ndat; k++) {
            tnew += waux[k] * SQUARE(daux[k] - edis[k]);
            smid += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        if (*ordinal) {
            (void)smacofSSFStressFList(edis, fval, gval, what, rpow, ndat);
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = fval[k];
            }
            (void)smacofSSFStressMonotone(ndat, ties, snew, iind, jind, blks,
                                          edis, dhat, wght);
            double sum = 0.0;
            for (int k = 0; k < Ndat; k++) {
                sum += wght[k] * SQUARE(dhat[k]);
            }
            sum = sqrt(sum);
            for (int k = 0; k < Ndat; k++) {
                dhat[k] /= sum;
            }
            *snew = 0.0;
            for (int k = 0; k < Ndat; k++) {
                *snew += wght[k] * SQUARE(dhat[k] - fval[k]);
            }
        } else {
            *snew = smid;
        }
        double ops = 0.0;
        for (int i = 0; i < Nobj * Ndim; i++) {
            ops = fmax(ops,
                       fmin(fabs(xold[i] - xnew[i]), fabs(xold[i] + xnew[i])));
        }
        if (*verbose) {
            if (*ordinal) {
                printf(
                    "itel %4d sold %*.*f told %*.*f tnew %*.*f smid %*.*f snew "
                    "%*.*f ops %*.*f\n",
                    *itel, *width, *digits, *sold, *width, *digits, told,
                    *width, *digits, tnew, *width, *digits, smid, *width,
                    *digits, *snew, *width, *digits, ops);
            } else {
                printf(
                    "itel %4d sold %*.*f told %*.*f tnew %*.*f  snew %*.*f ops "
                    "%*.*f\n",
                    *itel, *width, *digits, *sold, *width, *digits, told,
                    *width, *digits, tnew, *width, *digits, *snew, *width,
                    *digits, ops);
            }
        }
        if ((*itel == *itmax) || (ops < Eps)) {
            break;
        }
        for (int k = 0; k < Nobj * Ndim; k++) {
            xold[k] = xnew[k];
        }
        *sold = *snew;
        *itel += 1;
    }
    xfree(vinv);
    xfree(waux);
    xfree(daux);
    xfree(fval);
    xfree(gval);
    return;
}

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

void smacofSSFStressFList(double* x, double* f, double* g, int* what, double* rpow,
                          int* nvec) {
    double Rpow = *rpow;
    int Nvec = *nvec, What = *what;
    if (What == 0) {
        for (size_t i = 0; i < Nvec; i++) {
            f[i] = pow(x[i], Rpow);
            g[i] = Rpow * pow(x[i], Rpow - 1.0);
        }
        return;
    }
    if (What == 1) {
        for (size_t i = 0; i < Nvec; i++) {
            f[i] = log(x[i]);
            g[i] = 1.0 / x[i];
        }
        return;
    }
    if (What == 2) {
        for (size_t i = 0; i < Nvec; i++) {
            f[i] = g[i] = exp(x[i]);
        }
        return;
    }
}