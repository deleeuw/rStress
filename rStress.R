rStress <- function(delta,
                    xinit,
                    wght = 1 - diag(nrow(xinit)),
                    rpow = 1,
                    dpow = FALSE,
                    itmax = 10000,
                    eps = 1e-6,
                    verbose = TRUE) {
  itel <- 1
  nobj <- nrow(xinit)
  mobj <- 1 / nobj
  if (dpow) {
    dhat <- delta^rpow
  } else {
    dhat <- delta
  }
  xold <- xinit
  dold <- as.matrix(dist(xold))
  sold <- sum(wght * (dhat - (dold^rpow))^2) / 2
  repeat {
    wmat <- (rpow^2) * wght * ((dold + diag(nobj))^(2 * (rpow - 1)))
    dmat <- (dhat + (rpow - 1) * (dold^rpow)) / (rpow * (dold + diag(nobj))^(rpow - 1))
    hmat <- 0.5 * (dmat < 0) * (1 / (dold + diag(nobj)))
    sgno <- sum(wmat * (dmat - dold)^2) / 2 # always sgno = sold
    vmat <- -(wmat + hmat)
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + mobj) - mobj
    bmat <- -(dmat > 0)  * wmat * dmat / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% bmat %*% xold
    dnew <- as.matrix(dist(xnew))
    apse <- max(abs(dold - dnew))
    sgnn <- sum(wmat * (dmat - dnew)^2) / 2 # always sgnn < sgno = sold
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
        " sgno",
        formatC(
          sgno,
          digits = 8,
          width = 10,
          format = "f"
        ),
        " sgnn",
        formatC(
          sgnn,
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
  return(list(x = xnew, d = dnew, loss = snew, itel = itel))
}