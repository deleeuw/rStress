source("smacofDataUtilities.R")

smallMatrix <- matrix(c(
  0,1,2,4,7,11,
  1,0,3,5,8,12,
  2,3,0,6,9,13,
  4,5,6,0,10,14,
  7,8,9,10,0,15,
  11,12,13,14,15,0), 6, 6)
smallWeights <- matrix(c(
  0,1,2,3,4,5,
  1,0,2,3,4,5,
  2,2,0,3,4,5,
  3,3,3,0,4,5,
  4,4,4,4,0,5,
  5,5,5,5,5,0), 6, 6)
smallMatrix <- smallMatrix * sqrt(2 / sum(smallMatrix^2))
smallDist <- as.dist(smallMatrix)
smallData <- makeMDSData(smallDist, as.dist(smallWeights))
