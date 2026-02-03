
source("smacofDataUtilities.R")

data(ekman, package = "smacof")
ekmanLabels <- as.character(attr(ekman, "Labels"))
ekmanDist <- 1 - ekman
ekmanMatrix <- as.matrix(ekmanDist)
ekmanData <- makeMDSData(ekmanDist)
ikman <- ekmanMatrix^2
rkman <- apply(ikman, 1, mean)
mkman <- mean(ikman)
ckman <- -(ikman - outer(rkman, rkman, "+") + mkman) / 2
skman <- eigen(ckman)
xinit <- skman$vectors[, 1:2] %*% diag(sqrt(skman$values[1:2]))


