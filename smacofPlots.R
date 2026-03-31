source("smacofSSFStressSelect.R")

smacofShepardPlot <-
  function(h,
           fitlines = TRUE,
           colline = "RED",
           colpoint = "BLUE",
           lwd = 2,
           cex = 1,
           pch = 16) {
    what <- h$what
    rpow <- h$rpow
    hfnc <- smacofSSFStressSelect(what)
    del <- h$delt
    dis <- h$dist
    dht <- h$dhat
    x <- switch(h$usef + 1, del, hfnc$func(del, rpow))
    xlab <- switch(h$usef + 1, "delta", "f(delta)")
    y <- dht
    z <- hfnc$func(dis, rpow)
    ylab <- "dhat and f(dist)"
    main <- paste("ShepardPlot ", h$labl, ", what = ", what, ", rpow = ", rpow, sep = "")
    if (h$ordi) {
      main <- paste(main, ", ordinal", sep = "")
    } else {
      main <- paste(main, ", numerical", sep = "")
    }
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(min(x), max(x)),
      ylim = c(min(c(z, y)), max(c(z, y))),
      xlab = xlab,
      ylab = ylab,
      main = main,
      type = "n"
    )
    points(x,
           z,
           col = colpoint,
           cex = cex,
           pch = pch)
    points(x,
           y,
           col = colline,
           cex = cex,
           pch = pch)
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    lines(x,
          y,
          type = "l",
          lwd = lwd,
          col = colline)
  }

smacofConfigurationPlot <-
  function(h,
           labels = NULL,
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.0) {
    xmat <- h$xmat
    main <- paste("Configuration Plot ",
                  h$labl,
                  ", what = ",
                  what,
                  ", rpow = ",
                  rpow,
                  sep = "")
    if (h$ordi) {
      main <- paste(main, ", ordinal", sep = "")
    } else {
      main <- paste(main, ", numerical", sep = "")
    }
    if (is.null(labels)) {
      plot(
        xmat[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        pch = pch,
        col = col,
        cex = cex
      )
    }
    else {
      plot(
        xmat[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        type = "n"
      )
      text(xmat[, c(dim1, dim2)], labels, col = col, cex = cex)
    }
  }

smacofDistDhatPlot <- function(h,
                               fitlines = TRUE,
                               colline = "RED",
                               colpoint = "BLUE",
                               cex = 1.0,
                               lwd = 2,
                               pch = 16) {
  what <- h$what
  rpow <- h$rpow
  hfnc <- smacofSSFStressSelect(what)
  x <- hfnc$func(h$dist, rpow)
  y <- h$dhat
  xlab <- "f(dist)"
  ylab <- "dhat"
  main <- paste("DistDhatPlot ", h$labl, ", what = ", what, ", rpow = ", rpow, sep = "")
  if (h$ordi) {
    main <- paste(main, ", ordinal", sep = "")
  } else {
    main <- paste(main, ", numerical", sep = "")
  }
  uppe <- max(c(x, y))
  plot(
    x,
    y,
    xlab = xlab,
    ylab = ylab,
    xlim = c(0, uppe),
    ylim = c(0, uppe),
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = colline, lwd = lwd)
  if (fitlines) {
    m <- length(x)
    for (i in 1:m) {
      xx <- x[i]
      yy <- y[i]
      zz <- (xx + yy) / 2
      a <- matrix(c(xx, zz, yy, zz), 2, 2)
      lines(a, lwd = lwd)
    }
  }
}

smacofResidualPlot <- function(h,
                               probs = seq(0, 1, 0.25),
                               whatres = 1,
                               whatplot = 1,
                               qlines = TRUE,
                               colpoints = "RED",
                               collines = "BLUE",
                               lwdpoints = 2,
                               lwdlines = 2,
                               cex = 1,
                               pch = 16) {
  what <- h$what
  rpow <- h$rpow
  hfnc <- smacofSSFStressSelect(what)
  res <- hfnc$func(h$dist, rpow) - h$dhat
  res <- switch(whatres, res, abs(res), res^2)
  q <- quantile(res, probs)
  n <- length(q)
  e <- switch(whatplot, ecdf(res), density(res))
  main <- paste("ResidualPlot ", h$labl, ", what = ", what, ", rpow = ", rpow, sep = "")
  if (h$ordi) {
    main <- paste(main, ", ordinal", sep = "")
  } else {
    main <- paste(main, ", numerical", sep = "")
  }
  plot(e,
       col = colpoints,
       main = main,
       lwd = lwdpoints)
  for (i in 1:n) {
    abline(v = q[i], col = collines, lwd = lwdlines)
  }
}