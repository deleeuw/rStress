
source("smacofDataUtilities.R")

data(ekman, package = "smacof")
ekmanLabels <- as.character(attr(ekman, "Labels"))
ekman <- as.matrix(1 - ekman)
ikman <- ekman^2
rkman <- apply(ikman, 1, mean)
mkman <- mean(ikman)
ckman <- -(ikman - outer(rkman, rkman, "+") + mkman) / 2
skman <- eigen(ckman)
xini <- skman$vectors[, 1:2] %*% diag(sqrt(skman$values[1:2]))


