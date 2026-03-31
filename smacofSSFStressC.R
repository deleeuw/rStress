
dyn.load("smacofSSFStress.so")

suppressPackageStartupMessages(library(MASS, quietly = TRUE))
suppressPackageStartupMessages(library(isotone, quietly = TRUE))

source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

smacofSSFStressC <- function(theData,
                            ndim = 2,
                            xinit = NULL,
                            ties = 1,
                            itmax = 10000,
                            usef = 1,
                            eps = 1e-6,
                            what = 0,
                            rpow = 1,
                            digits = 8,
                            width = 10,
                            verbose = TRUE,
                            ordinal = FALSE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  blks <- theData$blocks
  delt <- theData$delta
  wght <- theData$weights / sum(theData$weights)
  labl <- theData$label
  if (is.null(xinit)) {
    xold <- smacofTorgerson(theData, ndim)$conf
  } else {
    xold <- xinit
  }
  dold <- dnew <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    dold[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  h <- smacofSSFStressSelect(what)
  func <- h$func
  if (usef == 1) {
    dhat <- func(delt, rpow)
  } else {
    dhat <- delt
  }
  if (ordinal) {
    dhat <- dhat / sqrt(sum(wght * dhat^2))
  }
  labd <- sum(wght * dold * dhat) / sum(wght * dold ^ 2)
  xold <- xold * labd
  dold <- dold * labd
  sold <- sum(wght * (dhat - func(dold, rpow))^2)
  itel <- 1
  snew <- 0.0
  xnew <- xold
  h <- .C(
    "smacofSSFStressEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    ties = as.integer(ties),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    ordinal = as.integer(ordinal),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    what = as.integer(what),
    rpow = as.double(rpow),
    iind = as.integer(iind - 1),
    jind = as.integer(jind - 1),
    blks = as.integer(blks),
    wght = as.double(wght),
    dhat = as.double(dhat),
    dold = as.double(dold),
    dnew = as.double(dnew),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  result <- list(
    delt = theData$delta,
    dhat = h$dhat,
    dist = h$dnew,
    xmat = matrix(h$xnew, nobj, ndim),
    wmat = h$wght,
    loss = h$snew,
    ndim = ndim,
    init = xinit,
    itel = h$itel,
    nobj = nobj,
    iind = h$iind,
    jind = h$jind,
    ordi = ordinal,
    ties = h$ties,
    what = what,
    rpow = rpow,
    labl = labl,
    usef = usef,
    date = date()
  )
  return(result)
}
