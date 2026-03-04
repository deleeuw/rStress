
library(RSpectra)

smacofTorgerson <- function(theData,  ndim = 2, itmax = 5, eps = 1e-6, verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delta <- theData$delta
  dmat <- matrix(0, nobj, nobj) # order n matrix of dissimilarities
  indi <- 1 - diag(nobj) # order n matrix with missing indicator
  mdel <- mean(theData$delta)
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- delta[k]
    indi[i, j] <- indi[j, i] <- 0
  }
  nmis <- sum(indi) / 2
  kind <- matrix(0, nmis, 2)
  k <- 1
  for (i in 2:nobj) {
    for (j in 1:(i - 1)) {
      if (indi[i,j] == 1) {
        kind[k, 1] <- i
        kind[k, 2] <- j
        k <- k + 1
      }
    }
  }
  cmat <- matrix(0, nmis, nmis)
  for (k in 1:nmis) {
    ij <- kind[k, ]
    for (l in 1:nmis) {
      kl <- kind[l, ]
      si <- sum(outer(ij, kl, "==")) + 1
      cmat[k, l] <- cmat[l, k] <- switch(si, 4, 4 - 2 * n, 4 * n^2 - 4 * n + 4) / n^2
    }
  }
  dmat <- dmat^2
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  bmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigs_sym(bmat, ndim, which = "LA")
  evev <- diag(sqrt(pmax(0, ev$values)))
  x <- ev$vectors[, 1:ndim] %*% evev
}
