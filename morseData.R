source("smacofDataUtilities.R")

data(morse, package = "smacof")
morseLabels <- as.character(attr(morse, "Labels"))
morseData <- makeMDSData(morse, label = "morseData")
morseMatrix <- as.matrix(morse)
imorse <- morseMatrix^2
rmorse <- apply(imorse, 1, mean)
mmorse <- mean(imorse)
cmorse <- -(imorse - outer(rmorse, rmorse, "+") + mmorse) / 2
smorse <- eigen(cmorse)
xini <- smorse$vectors[, 1:2] %*% diag(sqrt(smorse$values[1:2]))
