

library(MASS)
source("smacofDataUtilities.R")
source("smacofSSFStressR.R")

func <- function(x)
  log(x + 1)
gunc <- function(x)
  1 / (x + 1)
w0 <- as.dist(matrix(c(0, 2, 3, 2, 0, 4, 3, 4, 0), 3, 3))
d0 <- as.dist(matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3, 3))
x0 <- matrix(c(1, 0, -1, 0, 1, -1), 3, 2)

tester <- function(delt = d0,
                   wght = w0,
                   xinit = x0,
                   f = func,
                   g = gunc) {
  wght <- wght / sum(wght)
  xold <- xinit
  dold <- dist(xold)
  labd <- sum(wght * delt * dold) / sum(wght * dold^2)
  xold <- xold * labd
  dold <- dold * labd
  fdelt <- f(delt)
  dhat <- fdelt / sqrt(sum(wght * fdelt^2))
  fdold <- f(dold)
  sold <- sum(wght * (dhat - fdold)^2)
  gdold <- g(dold)
  waux <- wght * gdold^2
  daux <- (dhat - fdold) / gdold + dold
  told <- sum(waux * (daux - dold)^2)
  baux <- waux * daux / dold
  bmat <- -as.matrix(baux)
  diag(bmat) <- -rowSums(bmat)
  vmat <- -as.matrix(waux)
  diag(vmat) <- -rowSums(vmat)
  vinv <- ginv(vmat)
  xtmp <- bmat %*% xold
  xnew <- vinv %*% xtmp
  dnew <- dist(xnew)
  fdnew <- f(dnew)
  tnew <- sum(waux * (daux - dnew)^2)
  snew <- sum(wght * (dhat - fdnew)^2)
  print(c(sold, told, tnew, snew))
}
