source("smacofDataUtilities.R")

data(morse, package = "smacof")
morseLabels <- as.character(attr(morse, "Labels"))
morse <- as.matrix(morse)
imorse <- morse^2
rmorse <- apply(imorse, 1, mean)
mmorse <- mean(imorse)
cmorse <- -(imorse - outer(rmorse, rmorse, "+") + mmorse) / 2
smorse <- eigen(cmorse)
xini <- smorse$vectors[, 1:2] %*% diag(sqrt(smorse$values[1:2]))
