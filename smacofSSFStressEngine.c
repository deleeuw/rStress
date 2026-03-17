
#include "smacofSSFStress.h"

void smacofSSRStressEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal, double* sold,
                           double* snew, double* eps, int* what, double* rpow,
                           int* iind, int* jind, int* blks, double* wght,
                           double* dhat, double* dold, double* dnew,
                           double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double Eps = *eps;
    double* waux = xmalloc(Ndat * sizeof(double));
    double* daux = xmalloc(Ndat * sizeof(double));
    double* vaux = xmalloc(Ndat * sizeof(double));
    double* baux = xmalloc(Ndat * sizeof(double));
    double* fval = xmalloc(Ndat * sizeof(double));
    double* gval = xmalloc(Ndat * sizeof(double));
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    while (true) {
        double told = 0.0;
        (void)smacofSSFStressFList(dold, fval, gval, what, rpow, ndat);
        for (int k = 0; k < Ndat; k++) {
            waux[k] = wght[k] * SQUARE(gval[k]);
            daux[k] = (dhat[k] - fval[k]) / gval[k] + dold[k];
            if (daux[k] > 0) {
                vaux[k] = waux[k];
                baux[k] = waux[k] * daux[k] / dold[k];
            } else {
                vaux[k] = waux[k] - waux[k] * daux[k] / dold[k];
            }
            told += waux[k] * SQUARE(daux[k] - dold[k]);
        }
        (void)smacofMPInverseV(nobj, ndat, iind, jind, vaux, vinv);
        (void)smacofSSFStressMajorize(nobj, ndim, ndat, snew, iind, jind, baux,
                                      vinv, dhat, dold, dnew, xold, xnew);
        double tnew = 0.0, smid = 0.0;
        (void)smacofSSFStressFList(dnew, fval, gval, what, rpow, ndat);
        for (int k = 0; k < Ndat; k++) {
            tnew += waux[k] * SQUARE(daux[k] - dnew[k]);
            smid += wght[k] * SQUARE(dhat[k] - fval[k]);
        }
        if (*ordinal) {
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = fval[k];
            }
            (void)smacofSSFStressMonotone(ndat, ties, snew, iind, jind, blks,
                                          fval, dhat, wght);
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
        for (int k = 0; k < Ndat; k++) {
            ops = fmax(ops, fabs(dold[k] - dnew[k]));
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
        for (int k = 0; k < Ndat; k++) {
            dold[k] = dnew[k];
        }
        *sold = *snew;
        *itel += 1;
    }
    xfree(vinv);
    xfree(vaux);
    xfree(baux);
    xfree(daux);
    xfree(fval);
    xfree(gval);
    return;
}
