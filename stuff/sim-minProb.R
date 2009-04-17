
### If  D_n = D(X_1, ... X_n) is the dip statistic, we have
###     D_n >= 1/(2 n)
###   pm(n) := Pr[ D == 1/(2 n) ]  > 0
###   ==      ------------------------
## and hence the distribution of D_n(X), X ~ U[0,1] starts with a jump,
## from 0 to pm(n) at  d = 1/(2n).

### Now estimate pm(n) via simulation.
### The simulation is analyzed in file ./minProb-anal.R
###                                      ~~~~~~~~~~~~~~
setwd("~/R/Pkgs/diptest/stuff")

library(diptest)

Ns <- 500000 # number of samples (per n)
isim <- 1:Ns
nn <- c(4:18,20,25,30,35,40,50,60, 75, 100)

nMin <- sapply(nn, function(n)
               sum(sapply(isim,
                          function(i) abs(1 - 2*n*dip(runif(n))) < 1e-5)))
names(nMin) <- paste(nn)
attr(nMin, "Ns") <- Ns
## nMin / Ns == pm(n), i.e. pm(nn)
save(nMin, file= "minProb.rda")

proc.time()
sessionInfo()
