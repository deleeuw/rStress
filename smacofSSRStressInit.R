
library(RSpectra)

smacofSSRStressInit <- function(theData,  ndim = 2, rpow = 1) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  wght <- theData$weights
  dmat <- matrix(0, nobj, nobj)
  dhat <- (theData$delta)^2
  mdel <- mean(dhat)
  dmat <- mdel * (1 - diag(nobj))
  wmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- dhat[k]
    wmat[i, j] <- wmat[j, i] <- wght[k]
  }
  wsum <- sum(wmat)
  dmat <- dmat * sqrt(wsum / sum(wmat * dmat^2))
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  bmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigs_sym(bmat, ndim, which = "LA")
  evev <- diag(sqrt(pmax(0, ev$values)))
  x <- ev$vectors[, 1:ndim] %*% evev
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    edis[k] <- sqrt(sum((x[i,] - x[j, ])^2))
  }
  dhat <- theData$delta^rpow
  dhat <- dhat / sqrt(sum(wght * dhat ^ 2))
  dpow <- edis^rpow
  sdd <- sum(wght * dpow^2)
  sde <- sum(wght * dhat * dpow)
  lbd <- (sde / sdd)^(1 / rpow)
  edis <- lbd * edis
  x <- lbd * x
  rstress <- sum(wght * (dhat - edis^rpow)^2)
  dmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- edis[k]
  }
  return(list(x = x, edis = edis, dhat = dhat, rstress = rstress))
}

smacofSSRStressScale <- function(theData, xinit, rpow) {
  ndat <- theData$ndat
  wght <- theData$weights
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    edis[k] <- sqrt(sum((x[i,] - x[j, ])^2))
  }
  dhat <- theData$delta^rpow
  dhat <- dhat / sqrt(sum(wght * dhat ^ 2))
  dpow <- edis^rpow
  sdd <- sum(wght * dpow^2)
  sde <- sum(wght * dhat * dpow)
  lbd <- (sde / sdd)^(1 / rpow)
  edis <- lbd * edis
  x <- lbd * x
  rstress <- sum(wght * (dhat - edis^rpow)^2)
  result <- smacofSSRStressScale(theData)
  return(result)
}
