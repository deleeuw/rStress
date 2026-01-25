fStress <- function(delta,
                    xini,
                    wght = 1 - diag(nrow(xini)),
                    f = fPower1,
                    g = gPower1,
                    itmax = 10000,
                    eps = 1e-6,
                    verbose = TRUE) {
  itel <- 1
  nobj <- nrow(xini)
  mobj <- 1 / nobj
  delf <- f(delta)
  dini <- as.matrix(dist(xini))
  fini <- f(dini)
  labd <- sum(wght * delf * fini) / sum(wght * (fini^2))
  xold <- labd * xini
  dold <- as.matrix(dist(xold))
  sold <- sum(wght * (delf - f(dold))^2)
  repeat {
    wmat <- wght * g(dold)^2
    dmat <- ((delf - f(dold)) / g(dold)) + dold
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
    snew <- sum(wght * (delf - f(dnew))^2)
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
    x = xnew,
    d = dnew,
    loss = snew,
    itel = itel
  ))
}

fPower2 <- function(x) {
  nobj <- nrow(x)
  for (i in 1:nobj) {
    for (j in 1:nobj) {
      if (i == j) {
        next
      }
      x[i, j] <- x[i, j] ^ 2
    }
  }
  return(x)
}

gPower2 <- function(x) {
  nobj <- nrow(x)
  return((2 * x) + diag(nobj))
}

fPower1 <- function(x) {
  return(x)
}

gPower1 <- function(x) {
    nobj <- dim(x)[1]
    return(matrix(1, nobj, nobj))
}

fPowerh <- function(x) {
  return(sqrt(x))
}

gPowerh <- function(x) {
  nobj <- nrow(x)
  return(1 / sqrt(x + diag(nobj)) - diag(nobj))
}

flog <- function(x) {
  nobj <- nrow(x)
  return(log(x + diag(nobj)))
}

glog <- function(x) {
  nobj <- nrow(x)
  return(1 / (x + diag(nobj)))
}


fPower <- function(x, ...) {
  power <- function(x, r) {
    return(x ^ r)
  }
  r <- list(...)$r
  return(power(x, r))
}

gPower <- function(x, ...) {
  power <- function(x, r) {
    return(r * x ^ (r - 1))
  }
  r <- list(...)$r
  return(power(x, r))
}
  