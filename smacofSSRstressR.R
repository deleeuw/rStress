source("smacofDataUtilities.R")
source("smacofSSRStressInit.R")

smacofSSRStressR <- function(theData,
                             ndim = 2,
                             xinit = NULL,
                             rpow = 1,
                             itmax = 10000,
                             eps = 1e-6,
                             verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  h <- fromMDSData(theData)
  delta <- as.matrix(h$delta)
  wght <- as.matrix(h$weights)
  wght <- wght / sum(wght)
  dhat <- delta^rpow
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  if (is.null(xinit)) {
    hini <- smacofSSRStressInit(theData, ndim, rpow)
  } else {
    hini <- smacofSSRStressScale(theData, xinit, rpow)
  }
  xold <- matrix(hini$x, nobj, ndim)
  dold <- as.matrix(dist(xold))
  sold <- sum(wght * (dhat - dold^rpow)^2)
  mobj <- 1 / nobj
  itel <- 1
  repeat {
    wmat <- (rpow^2) * wght * ((dold + diag(nobj))^(2 * (rpow - 1)))
    dmat <- (dhat + (rpow - 1) * (dold^rpow)) / (rpow * (dold + diag(nobj))^(rpow - 1))
    hmat <- (dmat < 0) * (1 / (dold + diag(nobj)))
    told <- sum(wmat * (dmat - dold)^2) / 2 # always sgno = sold
    vmat <- -(wmat + hmat)
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + mobj) - mobj
    bmat <- -(dmat > 0)  * wmat * dmat / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% bmat %*% xold
    dnew <- as.matrix(dist(xnew))
    apse <- max(abs(dold - dnew))
    tnew <- sum(wmat * (dmat - dnew)^2) / 2 # always sgnn < sgno = sold
    snew <- sum(wght * (dhat - (dnew^rpow))^2) / 2
    if (verbose) {
      cat(
        " itel",
        formatC(itel, format = "d"),
        " sold",
        formatC(
          sold,
          digits = 8,
          width = 10,
          format = "f"
        ),
        " told",
        formatC(
          told,
          digits = 8,
          width = 10,
          format = "f"
        ),
        " tnew",
        formatC(
          tnew,
          digits = 8,
          width = 10,
          format = "f"
        ),
        " snew",
        formatC(
          snew,
          digits = 8,
          width = 10,
          format = "f"
        ),
        " apse",
        formatC(
          apse,
          digits = 8,
          width = 10,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (apse < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    dold <- dnew
    sold <- snew
  }
  return(list(
    x = xnew,
    d = dnew,
    loss = snew,
    itel = itel
  ))
}