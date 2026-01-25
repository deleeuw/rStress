rStress <- function(delta,
                    xini,
                    wght = 1 - diag(nrow(xini)),
                    r = 2,
                    itmax = 10000,
                    eps = 1e-6,
                    verbose = TRUE) {
  itel <- 1
  nobj <- nrow(xini)
  mobj <- 1 / nobj
  dpow <- delta^r
  dini <- as.matrix(dist(xini))
  eold <- dini^r
  labd <- (sum(wght * dpow * eold) / sum(wght * (eold^2)))^(1 / r)
  xold <- labd * xini
  dold <- as.matrix(dist(xold))
  sold <- sum(wght * (dpow - (dold^r))^2)
  repeat {
    wmat <- (r^2) * wght * ((dold + diag(nobj))^(2 * (r - 1)))
    dmat <- (dpow + (r - 1) * (dold^r)) / (r * (dold + diag(nobj))^(r - 1))
    sgno <- sum(wmat * (dmat - dold)^2) # always sgno = sold
    vmat <- -wmat
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + mobj) - mobj
    bmat <- -wmat * dmat / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% bmat %*% xold
    apse <- max(abs(xold - xnew))
    dnew <- as.matrix(dist(xnew))
    sgnn <- sum(wmat * (dmat - dnew)^2) # always sgnn < sgno = sold
    snew <- sum(wght * (dpow - (dnew^r))^2)
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