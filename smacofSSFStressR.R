

source("smacofTorgerson.R")
source("smacofSSFStressSelect.R")

suppressPackageStartupMessages(library(MASS, quietly = TRUE))
suppressPackageStartupMessages(library(isotone, quietly = TRUE))

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
                             ordinal = FALSE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  blks <- theData$blocks
  delt <- theData$delta
  labl <- theData$label
  wght <- theData$weights / sum(theData$weights)
  if (is.null(xinit)) {
    xold <- smacofTorgerson(theData, ndim)$conf
  } else {
    xold <- xinit
  }
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    dold[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  h <- smacofSSFStressSelect(what)
  func <- h$func
  gunc <- h$gunc
  dhat <- func(delt, rpow)
  if (ordinal) {
    dhat <- dhat / sqrt(sum(wght * dhat^2))
  }
  waux <- wght * gunc(delt, rpow)^2
  waux <- waux / sum(waux)
  labd <- sum(waux * dold * delt) / sum(waux * dold^2)
  xold <- xold * labd
  dold <- dold * labd
  sold <- sum(wght * (dhat - func(dold, rpow))^2)
  itel <- 1
  repeat {
    fval <- func(dold, rpow)
    gval <- gunc(dold, rpow)
    waux <- wght * gval^2
    daux <- ((dhat - fval) / gval) + dold
    told <- sum(waux * (daux - dold)^2)
    vmat <- bmat <- matrix(0, nobj, nobj)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      if (daux[k] > 0) {
        vmat[i, j] <- vmat[j, i] <- waux[k]
        bmat[i, j] <- bmat[j, i] <- waux[k] * daux[k] / dold[k]
      } else {
        vmat[i, j] <- vmat[j, i] <- waux[k] - (waux[k] * daux[k] / dold[k])
      }
    }
    vmat <- -vmat
    bmat <- -bmat
    diag(vmat) <- -rowSums(vmat)
    diag(bmat) <- -rowSums(bmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    xnew <- vinv %*% bmat %*% xold
    dnew <- rep(0, ndat)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      dnew[k] <- sqrt(sum((xnew[i, ] - xnew[j, ])^2))
    }
    ops <- max(abs(dold - dnew))
    fdis <- func(dnew, rpow)
    smid <- sum(wght * (dhat - fdis)^2)
    tnew <- sum(waux * (daux - dnew)^2)
    if (ordinal) {
      tchr <- switch(ties, "primary", "secondary", "tertiary")
      dhat <- gpava(
        z = delt,
        y = fdis,
        weights = wght,
        ties = tchr
      )$x
      dhat <- dhat / sqrt(sum(wght * dhat^2))
      snew <- sum(wght * (dhat - fdis)^2)
    } else {
      snew <- smid
    }
    if (verbose) {
      if (ordinal) {
        cat(
          "itel ",
          formatC(itel, digits = 4, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = digits,
            width = width,
            format = "f"
          ),
          "told ",
          formatC(
            told,
            digits = digits,
            width = width,
            format = "f"
          ),
          "tnew ",
          formatC(
            tnew,
            digits = digits,
            width = width,
            format = "f"
          ),
          "smid ",
          formatC(
            smid,
            digits = digits,
            width = width,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = digits,
            width = width,
            format = "f"
          ),
          "ops ",
          formatC(
            ops,
            digits = digits,
            width = width,
            format = "f"
          ),
          "\n"
        )
      } else {
        cat(
          "itel ",
          formatC(itel, digits = 4, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = digits,
            width = width,
            format = "f"
          ),
          "told ",
          formatC(
            told,
            digits = digits,
            width = width,
            format = "f"
          ),
          "tnew ",
          formatC(
            tnew,
            digits = digits,
            width = width,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = digits,
            width = width,
            format = "f"
          ),
          "ops ",
          formatC(
            ops,
            digits = digits,
            width = width,
            format = "f"
          ),
          "\n"
        )
      }
    }
    if ((itel == itmax) || (ops < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    itel <- itel + 1
  }
  result <- list(
    delt = delt,
    dhat = dhat,
    dist = dnew,
    xmat = xnew,
    wmat = wght,
    loss = snew,
    ndim = ndim,
    init = xinit,
    itel = itel,
    nobj = nobj,
    iind = iind,
    jind = jind,
    ordi = ordinal,
    ties = ties,
    what = what,
    rpow = rpow,
    labl = labl,
    date = date()
  )
  return(result)
}
