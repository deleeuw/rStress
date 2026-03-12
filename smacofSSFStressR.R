
source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

library(MASS)
library(isotone)

smacofSSFStressR <- function(theData,
                            ndim = 2,
                            xinit = NULL,
                            ties = 1,
                            itmax = 10000,
                            eps = 1e-6,
                            what = 0,
                            rpow = 1,
                            digits = 8,
                            width = 10,
                            verbose = TRUE,
                            weighted = FALSE,
                            ordinal = FALSE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  blks <- theData$blocks
  delt <- theData$delta
  wght <- theData$weights / sum(theData$weights)
  if (what == 0) {
    func <- function(x, r) x ^ r
    gunc <- function(x, r) r * x ^ (r - 1)
  }
  if (what == 1) {
    func <- function(x, r) log(r * x + 1)
    gunc <- function(x, r) r / log(r * x + 1)
  }
  if (what == 2) {
    func <- function(x, r) exp(r * x) - 1
    gunc <- function(x, r) r * exp(r * x)
  }
  if (is.null(xinit)) {
    xold <- smacofTorgerson(theData, ndim)$conf
  } 
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    dold[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  labd <- sum(wght * dold * delt) / sum(dold ^ 2)
  xold <- xold * labd
  dold <- dold * labd
  dhat <- func(delt, rpow)
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  sold <- sum(wght * (dhat - func(dold, rpow))^2)
  itel <- 1
  repeat {
    waux <- wght * gunc(dold, rpow)^2
    daux <- (dhat - func(dold, rpow)) / gunc(dold, rpow) + dold
    told <- sum(waux * (daux - dold)^2)
    vmat <- bmat <- matrix(0, nobj, nobj)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      if (daux[k] > 0) {
        vmat[i, j] <- vmat[j, i] <- waux[k]
        bmat[i, j] <- bmat[j, i] <- waux[k] * daux[k] / dold[k]
      } else {
        vmat[i, j] <- vmat[j, i] <- waux[k] - waux[k] * daux[k] / dold[k]
      }
    }
    vmat <- -vmat
    bmat <- -bmat
    diag(vmat) <- -rowSums(vmat)
    diag(bmat) <- -rowSums(bmat)
    vinv <- ginv(vmat)
    xnew <- vinv %*% bmat %*% xold
    epse <- max(abs(xold - xnew))
    dnew <- rep(0, ndat)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      dnew[k] <- sqrt(sum((xnew[i, ] - xnew[j, ])^2))
    }
    fdis <- func(dnew, rpow)
    smid <- sum(wght * (dhat - fdis)^2)
    tnew <- sum(waux * (daux - dnew)^2)
    if (ordinal) {
      tchr <- switch(ties, "primary", "secondary", "tertiary")
      dhat <- gpava(z = delt, y = fdis, weights = wght, ties = tchr)$x
      dhat <- dhat / sqrt(sum(dhat^2))
      snew <- sum(wght * (dhat - fdis)^2)
    } else {
      snew <- smid
    }
    if (verbose) {
      cat("itel ", formatC(itel, digits = 4, format = "d"),
          "sold ", formatC(sold, digits = digits, width = width, format = "f"),
          "told ", formatC(told, digits = digits, width = width, format = "f"),
          "tnew ", formatC(tnew, digits = digits, width = width, format = "f"),
          "snew ", formatC(snew, digits = digits, width = width, format = "f"),
          "epse ", formatC(epse, digits = digits, width = width, format = "f"),
          "\n")
    }
    if ((itel == itmax) || (epse < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    itel <- itel + 1
  }
  result <- list(
    delta = delt,
    dhat = dhat,
    confdist = dnew,
    conf = xnew,
    weightmat = wght,
    stress = snew,
    ndim = ndim,
    init = xinit,
    niter = itel,
    nobj = nobj,
    iind = iind,
    jind = jind,
    weighted = weighted,
    ordinal = ordinal,
    ties = ties, 
    what = what,
    rpow = rpow
  )
  return(result)
}
