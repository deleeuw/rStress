dyn.load("smacofSSRStress.so")

source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofSSRStressInit.R")

smacofSSRStress <- function(theData,
                            ndim = 2,
                            xinit = NULL,
                            ties = 1,
                            itmax = 1000,
                            eps = 1e-6,
                            rpow = 1.0,
                            digits = 8,
                            width = 10,
                            verbose = TRUE,
                            weighted = FALSE,
                            ordinal = FALSE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  wght <- theData$weights
  dhat <- theData$delta^rpow
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  if (is.null(xinit)) {
    h <- smacofSSRStressInit(theData, ndim, rpow)
    xold <- h$x
    edis <- h$edis
    sold <- h$rstress
  }
  itel <- 1
  blks <- theData$blocks
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
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}
