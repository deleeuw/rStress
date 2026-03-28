source("smacofDataUtilities.R")

data(wish, package = "smacof")
wishDist <- 7.0 - wish
wishMatrix <- as.matrix(wishDist)
wishData <- makeMDSData(wishDist, label = "wishData")
wishLabels <- as.character(attr(wish, "Labels"))


