smacofSSFStressSelect <- function(what = 0) {
  if (what == 0) {
    func <- function(x, r)
      x^r
    gunc <- function(x, r)
      r * x^(r - 1.0)
    finv <- function(y, r)
      y^(1.0 / r)
  }
  if (what == 1) {
    func <- function(x, r)
      log(r * x + 1.0)
    gunc <- function(x, r)
      r / (r * x + 1.0)
    finv <- function(y, r)
      (exp(y) - 1.0) / r
  }
  if (what == 2) {
    func <- function(x, r)
      exp(r * x) - 1
    gunc <- function(x, r)
      r * exp(r * x)
    finv <- function(y, r)
      log(y + 1) / r
  }
  return(list(
    func = func,
    gunc = gunc,
    finv = finv
  ))
}