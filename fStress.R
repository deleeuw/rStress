source("smacofTorgerson.R")
source("ekmanData.R")
dyn.load("smacofSSRStressMonotone.so")

fStress <- function(theData,
                    xinit = NULL,
                    ndim = 2,
                    f = fPower,
                    g = gPower,
                    rpow = 1,
                    ties = 1,
                    ordinal = FALSE,
                    itmax = 10000,
                    eps = 1e-6,
                    verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  blks <- theData$blocks
  wght <- theData$weights / sum(theData$weights)
  dhat <- f(theData$delta, rpow)
  self <- sum(wght * (dhat^2))
  dhat <- dhat / sqrt(self)
  if (is.null(xinit)) {
    xold <- smacofTorgerson(theData, ndim)$conf
  } else {
    xold <- xinit
  }
  self <- sum(wght * (dhat^2))
  dhat <- dhat / sqrt(self)
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    dold[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  fold <- f(dold, rpow)
  sold <- sum(wght * (dhat - fold)^2)
  itel <- 1
  repeat {
    wmat <- bmat <- matrix(0, nobj, nobj)
    daux <- ((dhat - fold) / g(dold, rpow)) + dold
    waux <- wght * g(dold, rpow)^2
    told <- sum(waux * (daux - dold)^2) # always told = sold
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      wmat[i, j] <- wmat[j, i] <- waux[k]
      bmat[i, j] <- bmat[j, i] <- -waux[k] * daux[k] / dold[k]
    }
    vmat <- -wmat
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% bmat %*% xold
    apse <- max(abs(xold - xnew))
    dnew <- rep(0, ndat)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      dnew[k] <- sqrt(sum((xnew[i, ] - xnew[j, ])^2))
    }
    fnew <- f(dnew, rpow)
    tnew <- sum(waux * (daux - dnew)^2) # always tnew < told = sold
    smid <- sum(wght * (dhat - fnew)^2)
    if (ordinal) {
      dhat <- f(dnew, rpow)
      h <- .C(
        "smacofSSRStressMonotone",
        ndat = as.integer(ndat),
        ties = as.integer(ties),
        snew = as.double(0.0),
        iind = as.integer(iind),
        jind = as.integer(jind),
        blks = as.integer(blks),
        dnew = as.double(dnew),
        fnew = as.double(fnew),
        wght = as.double(wght)
      )
      dhat <- h$fnew
      self <- sum(wght * (dhat^2))
      dhat <- dhat / sqrt(self)
      snew <- sum(wght * (dhat - fnew)^2)
    } else {
      snew <- smid
    }
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
        " smid",
        formatC(
          smid,
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
    fold <- fnew
    sold <- snew
  }
  return(
    list(
      delta = theData$delta,
      dhat = dhat,
      weightmat = wght,
      conf = xnew,
      confdist = dnew,
      stress = snew,
      niter = itel
    )
  )
}

flog <- function(x, r) {
  return(log(x))
}

glog <- function(x, r) {
  return(1 / x)
}

fPower <- function(x, r) {
  return(x^r)
}

gPower <- function(x, r) {
  return(r * x^(r - 1))
}

fexp <- function(x, r) {
  return(exp(x))
}

gexp <- function(x, r) {
  return(exp(x))
}