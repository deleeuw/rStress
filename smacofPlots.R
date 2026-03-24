source("smacofSSFStressSelect.R")

smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           fitlines = TRUE,
           colline = "RED",
           colpoint = "BLUE",
           type = 0,
           lwd = 2,
           cex = 1,
           pch = 16) {
    what <- h$what
    rpow <- h$rpow
    hfnc <- smacofSSFStressSelect(what)
    del <- h$delta
    dis <- h$confdist
    dht <- h$dhat
    if (type == 0) {
      x <- hfnc$func(del, rpow)
      y <- dht
      z <- hfnc$func(dis, rpow)
      xlab <- "f(delta)"
      ylab <- "dhat and f(dist)"
    } else {
      x <- del
      y <- hfnc$finv(dht, rpow)
      z <- dis
      xlab <- "delta"
      ylab <- "finv(dhat) and dist"
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
           main = "ConfigurationPlot",
           labels = NULL,
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.0) {
    xnew <- h$conf
    if (is.null(labels)) {
      plot(
        xnew[, c(dim1, dim2)],
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
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        type = "n"
      )
      text(xnew[, c(dim1, dim2)], labels, col = col, cex = cex)
    }
  }

smacofDistDhatPlot <- function(h,
                               fitlines = TRUE,
                               colline = "RED",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               type = 0,
                               cex = 1.0,
                               lwd = 2,
                               pch = 16) {
  what <- h$what
  rpow <- h$rpow
  hfnc <- smacofSSFStressSelect(what)
  if (type == 0) {
    x <- hfnc$func(h$confdis, rpow)
    y <- h$dhat
    xlab <- "f(dist)"
    ylab <- "dhat"
  } else {
    x <- h$confdis
    y <- hfnc$finv(h$dhat, rpow)
    xlab <- "dist"
    ylab <- "finv(dhat)"
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
                               main = "ResidualPlot",
                               probs = seq(0, 1, 0.25),
                               dats = 1,
                               whatplot = "ecdf",
                               type = 0,
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
  if (type == 0) {
    res <- h$confdist - hfnc$finv(h$dhat, rpow)
  } else {
    res <- hfnc$func(h$confdist, rpow) - h$dhat
  }
  res <- switch(dats, res, abs(res), res^2)
  q <- quantile(res, probs)
  n <- length(q)
  e <- switch(whatplot, "ecdf" = ecdf(res), "density" = density(res))
  plot(e,
       col = colpoints,
       main = main,
       lwd = lwdpoints)
  for (i in 1:n) {
    abline(v = q[i], col = collines, lwd = lwdlines)
  }
}