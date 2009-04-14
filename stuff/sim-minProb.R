
### If  D = D(X_1, ... X_n) is the dip statistic, we have
###       D >= 1/(2 n)
###   pm(n) := Pr[ D == 1/(2 n) ]  > 0
###   ==      ------------------------
## and hence the distribution of D(X), X ~ U[0,1] starts with a jump,
## from 0 to pm(n) at  d = 1/(2n).

### Now estimate pm(n) via simulation :

### This simulation is analyzed in file ./minProb-anal.R
##                                        ~~~~~~~~~~~~~~

library(diptest)

Ns <- 100000 # number of samples (per n)
isim <- 1:Ns
nn <- c(4:20,25,30,40,50)

nMin <- sapply(nn, function(n)
               sum(sapply(isim,
                          function(i) abs(1 - 2*n*dip(runif(n))) < 1e-5)))
names(nMin) <- paste(nn)
attr(nMin, "Ns") <- Ns
## nMin / Ns == pm(n), i.e. pm(nn)
save(nMin, file= "minProb.rda")

q('no')
