fStress <- function(delta,
                    xini,
                    wght = 1 - diag(nrow(xini)),
                    f = fPower,
                    g = gPower,
                    r = 1,
                    itmax = 10000,
                    eps = 1e-6,
                    verbose = TRUE) {
  itel <- 1
  nobj <- nrow(xini)
  mobj <- 1 / nobj
  delf <- f(delta, r)
  self <- sum(wght * (delf^2))
  dini <- as.matrix(dist(xini))
  fini <- f(dini, r)
  labd <- sum(wght * delf * fini) / sum(wght * (fini^2))
  xold <- labd * xini
  dold <- as.matrix(dist(xold))
  sold <- sum(wght * (delf - f(dold, r))^2) / self
  repeat {
    wmat <- wght * g(dold, r)^2
    dmat <- ((delf - f(dold, r)) / (g(dold, r) + diag(nobj))) + dold
    sgno <- sum(wmat * (dmat - dold)^2) / self# always sgno = sold
    vmat <- -wmat
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + mobj) - mobj
    bmat <- -wmat * dmat / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% bmat %*% xold
    apse <- max(abs(xold - xnew))
    dnew <- as.matrix(dist(xnew))
    sgnn <- sum(wmat * (dmat - dnew)^2) / self # always sgnn < sgno = sold
    snew <- sum(wght * (delf - f(dnew, r))^2) / self
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
  return(list(
    delta = delta, 
    dhat = delf, 
    weightmat = wght,
    conf = xnew,
    confdist = dnew,
    stress = snew,
    niter = itel
  ))
}

flog <- function(x, r) {
  nobj <- nrow(x)
  return(log(x + diag(nobj)))
}

glog <- function(x, r) {
  nobj <- nrow(x)
  return(1 / (x + diag(nobj)))
}


fPower <- function(x, r) {
  return(x^r)
}

gPower <- function(x, r) {
  diag(x) <- 1
  aux <- r * x^(r - 1)
  diag(aux) <- 0
  return(aux)
}

fexp <- function(x, r) {
  aux <- exp(x)
  diag(aux) <- 0
  return(aux)
}

gexp <- function(x, r) {
  aux <- exp(x)
  diag(aux) <- 0
  return(aux)
}