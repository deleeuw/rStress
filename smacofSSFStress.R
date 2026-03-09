dyn.load("smacofSSRStress.so")

source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

smacofSSRStress <- function(theData,
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
  wght <- theData$weights / sum(theData$weights)
  if (what == 0) {
    func <- function(x, r) x ^ r
  }
  if (what == 1) {
    func <- function(x, r) log(x)
  }
  if (what == 2) {
    func <- function(x, r) exp(x)
  }
  dhat <- func(theData$delta, rpow)
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  if (is.null(xinit)) {
    xold <- smacofTorgerson(theData, ndim)$conf
  } 
  sold <- 0.0
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
    sold <- sold + wght[k] * (dhat[k] - func(edis[k], rpow))^2
  }
  # scale
  itel <- 1
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- xold
  h <- .C(
    "smacofSSRStressEngine",
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
    weighted = as.integer(weighted),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    what = as.integer(what),
    rpow = as.double(rpow),
    iind = as.integer(iind - 1),
    jind = as.integer(jind - 1),
    blks = as.integer(blks),
    wght = as.double(wght),
    edis = as.double(edis),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  result <- list(
    delta = theData$delta,
    dhat = h$dhat,
    confdist = h$edis,
    conf = matrix(h$xnew, nobj, ndim),
    weightmat = h$wght,
    stress = h$snew,
    ndim = ndim,
    init = xinit,
    niter = h$itel,
    nobj = nobj,
    iind = h$iind,
    jind = h$jind,
    weighted = weighted,
    ordinal = ordinal,
    ties = h$ties
  )
  return(result)
}
