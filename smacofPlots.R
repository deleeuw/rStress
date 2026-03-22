source("smacofSSFStressSelect.R")

smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           fitlines = TRUE,
           colline = "RED",
           colpoint = "BLUE",
           resolution = 100,
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
    if (type == 1) {
      x <- del
      y <- hfnc$finv(dht, rpow)
      z <- dis
      xlab <- "delta"
      ylab <- "finv(dhat) and dist"
    } else {
      x <- hfnc$func(del, rpow)
      y <- dht
      z <- hfnc$func(dis, rpow)
      xlab <- "f(delta)"
      ylab <- "dhat and f(dist)"
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
                               cex = 1.0,
                               lwd = 2,
                               pch = 16) {
  what <- h$what
  rpow <- h$rpow
  if (what == 0) {
    func <- function(x, r) x ^ r
  }
  if (what == 1) {
    func <- function(x, r) log(r * x + 1)
  }
  if (what == 2) {
    func <- function(x, r) exp(r * x) - 1
  }
  fdis <- func(h$confdist, rpow)
  uppe <- max(c(fdis, h$dhat))
  plot(
    fdis,
    h$dhat,
    xlab = "distance",
    ylab = "disparity",
    xlim = c(0, uppe),
    ylim = c(0, uppe),
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = colline, lwd = lwd)
  if (fitlines) {
    m <- length(fdis)
    for (i in 1:m) {
      x <- fdis[i]
      y <- h$dhat[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, lwd = lwd)
    }
  }
}

smacofResidualPlot <- function(h, 
                               main = "ResidualPlot",
                               probs = seq(0, 1, 0.25),
                               dats = 1,
                               type = "ecdf",
                               qlines = TRUE,
                               colpoints = "RED",
                               collines = "BLUE",
                               lwdpoints = 2,
                               lwdlines = 2,
                               cex = 1,
                               pch = 16) {
  res <- h$confdist - h$dhat
  res <- switch(dats, res, abs(res), res^2)
  q <- quantile(res, probs)
  n <- length(q)
  e <- switch(type, 
         "ecdf" = ecdf(res),
         "density" = density(res)
         )
  plot(e, col = colpoints, main = main, lwd = lwdpoints)
  for (i in 1:n) {
    abline(v = q[i], col = collines, lwd = lwdlines)
  }
}