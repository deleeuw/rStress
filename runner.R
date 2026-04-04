library(knitr)
library(microbenchmark)

source("smacofSSFStressC.R")
source("smacofSSFStressR.R")
source("smacofTorgerson.R")
source("smacofSSFStressSelect.R")
source("ekmanData.R")

dyn.load("smacofSSFStress.so")

hl <- smacofSSFStressR(ekmanData, what = 1, verbose = FALSE)
hp1 <- smacofSSFStressR(ekmanData, rpow = .1, verbose = FALSE)
hp5 <- smacofSSFStressR(ekmanData, rpow = .5, verbose = FALSE)
h1 <- smacofSSFStressR(ekmanData, rpow = 1, verbose = FALSE)
h2 <- smacofSSFStressR(ekmanData, rpow = 2, verbose = FALSE)
h5 <- smacofSSFStressR(ekmanData, rpow = 5, verbose = FALSE)
he <- smacofSSFStressR(ekmanData, what = 2, verbose = FALSE)

rloss <- c(hl$loss, hp1$loss, hp5$loss, h1$loss,h2$loss,h5$loss,he$loss)
riter <- c(hl$itel, hp1$itel, hp5$itel, h1$itel,h2$itel,h5$itel,he$itel)
x <- matrix(c(riter, rloss), 7, 2)
rnames <- c("logarithm", 
            "power .1",
            "power .5",
            "power 1",
            "power 2", 
            "power 5",
            "exponent"
)
row.names(x) <- rnames
cnames <- c("itel", "fStress")
colnames(x) <- cnames
print(x)

hlc <- smacofSSFStressC(ekmanData, what = 1, verbose = FALSE)
hp1c <- smacofSSFStressC(ekmanData, rpow = .1, verbose = FALSE)
hp5c <- smacofSSFStressC(ekmanData, rpow = .5, verbose = FALSE)
h1c <- smacofSSFStressC(ekmanData, rpow = 1, verbose = FALSE)
h2c <- smacofSSFStressC(ekmanData, rpow = 2, verbose = FALSE)
h5c <- smacofSSFStressC(ekmanData, rpow = 5, verbose = FALSE)
hec <- smacofSSFStressC(ekmanData, what = 2, verbose = FALSE)

closs <- c(hlc$loss, hp1c$loss, hp5c$loss, h1c$loss,h2c$loss,h5c$loss,hec$loss)
citer <- c(hlc$itel, hp1c$itel, hp5c$itel, h1c$itel,h2c$itel,h5c$itel,hec$itel)
x <- matrix(c(citer, closs), 7, 2)
rnames <- c("logarithm", 
            "power .1",
            "power .5",
            "power 1",
            "power 2", 
            "power 5",
            "exponent"
)
row.names(x) <- rnames
cnames <- c("itel", "fStress")
colnames(x) <- cnames
print(x)

hlo <- smacofSSFStressC(ekmanData, what = 1, verbose = FALSE, ordinal = TRUE)
hp1o <- smacofSSFStressC(ekmanData, rpow = .1, verbose = FALSE, ordinal = TRUE)
hp5o <- smacofSSFStressC(ekmanData, rpow = .5, verbose = FALSE, ordinal = TRUE)
h1o <- smacofSSFStressC(ekmanData, rpow = 1, verbose = FALSE, ordinal = TRUE)
h2o <- smacofSSFStressC(ekmanData, rpow = 2, verbose = FALSE, ordinal = TRUE)
h5o <- smacofSSFStressC(ekmanData, rpow = 5, verbose = FALSE, ordinal = TRUE)
heo <- smacofSSFStressC(ekmanData, what = 2, verbose = FALSE, ordinal = TRUE)

oloss <- c(hlo$loss, hp1o$loss, hp5o$loss, h1o$loss,h2o$loss,h5o$loss,heo$loss)
oiter <- c(hlo$itel, hp1o$itel, hp5o$itel, h1o$itel,h2o$itel,h5o$itel,heo$itel)
x <- matrix(c(oiter, oloss), 7, 2)
rnames <- c(
  "logarithm",
  "power .1",
  "power .5",
  "power 1",
  "power 2", 
  "power 5",
  "exponential"
)
row.names(x) <- rnames
cnames <- c("itel", "fStress")
colnames(x) <- cnames
print(x)

hmetric <- microbenchmark(smacofSSFStressR(ekmanData, what = 1, verbose = FALSE),
smacofSSFStressC(ekmanData, what = 1, verbose = FALSE),
smacofSSFStressR(ekmanData, rpow = 0.1, verbose = FALSE),
smacofSSFStressC(ekmanData, rpow = 0.1, verbose = FALSE),
smacofSSFStressR(ekmanData, rpow = 0.5, verbose = FALSE),
smacofSSFStressC(ekmanData, rpow = 0.5, verbose = FALSE),
smacofSSFStressR(ekmanData, rpow = 1, verbose = FALSE),
smacofSSFStressC(ekmanData, rpow = 1, verbose = FALSE),
smacofSSFStressR(ekmanData, rpow = 2, verbose = FALSE),
smacofSSFStressC(ekmanData, rpow = 2, verbose = FALSE),
smacofSSFStressR(ekmanData, rpow = 5, verbose = FALSE),
smacofSSFStressC(ekmanData, rpow = 5, verbose = FALSE),
smacofSSFStressR(ekmanData, what = 2, verbose = FALSE),
smacofSSFStressC(ekmanData, what = 2, verbose = FALSE))
print(hmetric)

hnmetric <- microbenchmark(smacofSSFStressR(ekmanData, what = 1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, what = 1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, rpow = 0.1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, rpow = 0.1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, rpow = 0.5, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, rpow = 0.5, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, rpow = 1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, rpow = 1, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, rpow = 2, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, rpow = 2, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, rpow = 5, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, rpow = 5, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressR(ekmanData, what = 2, verbose = FALSE, ordinal = TRUE),
                          smacofSSFStressC(ekmanData, what = 2, verbose = FALSE, ordinal = TRUE))
print(hnmetric)

dyn.unload("smacofSSFStress.so")
