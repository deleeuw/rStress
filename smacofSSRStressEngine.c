#include "smacofSSRStress.h"

void smacofSSRStressEngine(int* nobj, int* ndim, int* ndat, int* itel,
                           int* ties, int* itmax, int* digits, int* width,
                           int* verbose, int* ordinal, int* weighted,
                           double* sold, double* snew, double* eps,
                           double* rpow, int* iind, int* jind, int* blks,
                           double* wght, double* edis, double* dhat,
                           double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double Rpow = *rpow, Eps = *eps;
    double* wadj = xmalloc(Ndat * sizeof(double));
    double* dadj = xmalloc(Ndat * sizeof(double));
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    while (true) {
        double told = 0.0;
        for (int k = 0; k < Ndat; k++) {
            wadj[k] = SQUARE(Rpow) * wght[k] * pow(edis[k], 2.0 * (Rpow - 1.0));
            dadj[k] = (dhat[k] + (Rpow - 1.0) * pow(edis[k], Rpow)) /
                      (Rpow * pow(edis[k], Rpow - 1.0));
            told += wadj[k] * SQUARE(dadj[k] - edis[k]);
        }
        (void)smacofMPInverseV(nobj, ndat, iind, jind, wadj, vinv);
        (void)smacofSSRStressMajorize(nobj, ndim, ndat, snew, iind, jind,
                                      weighted, wadj, vinv, edis, dadj, xold,
                                      xnew);
        double tnew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            tnew += wadj[k] * SQUARE(dadj[k] - edis[k]);
        }
        double smid = smacofSSRStressLoss(ndat, edis, dhat, wght, rpow);
        if (*ordinal) {
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = pow(edis[k], Rpow);
            }
            (void)smacofSSRStressMonotone(ndat, ties, snew, iind, jind, blks,
                                          edis, dhat, wght);
            double sum = 0.0;
            for (int k = 0; k < Ndat; k++) {
                sum += wght[k] * SQUARE(dhat[k]);
            }
            sum = sqrt(sum);
            for (int k = 0; k < Ndat; k++) {
                dhat[k] /= sum;
            }
            *snew = smacofSSRStressLoss(ndat, edis, dhat, wght, rpow);
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
    xfree(wadj);
    xfree(dadj);
    return;
}

double smacofSSRStressLoss(int* ndat, double* edis, double* dhat, double* wght,
                           double* rpow) {
    double loss = 0.0;
    for (int k = 0; k < *ndat; k++) {
        loss += wght[k] * SQUARE(dhat[k] - pow(edis[k], *rpow));
    }
    return loss;
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