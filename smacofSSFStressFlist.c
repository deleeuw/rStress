#include "smacofSSFStress.h"

void smacofSSFStressFList(double* x, double* f, double* g, int* what,
                          double* rpow, int* nvec) {
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
            f[i] = log(Rpow * x[i] + 1.0);
            g[i] = Rpow / (Rpow * x[i] + 1.0);
        }
        return;
    }
    if (What == 2) {
        for (size_t i = 0; i < Nvec; i++) {
            f[i] = exp(Rpow * x[i]) - 1;
            g[i] = Rpow * exp(Rpow * x[i]);
        }
        return;
    }
    if (What == 3) {
      for (size_t i = 0; i < Nvec; i++) {
        f[i] = Rpow * log(x[i]);
        g[i] = Rpow / x[i];
      }
      return;
    }
}