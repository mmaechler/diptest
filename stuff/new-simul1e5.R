#### More extensive simulations than  Hartigan (1985)'s  table 1
#### ---------------------------------------------------
##-> ./dip-simul.S
##-> ./sim1.R
library(diptest)

nn <- c(4:10, 15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000)
Ns <- 100001
##    ~~~~~~ number of "simulations" (i.e. samples for each n in nn)

## which percentages
p.hi <- sort(1 - c(outer(c(1,2,5), 2:5, function(x,y) x*10^-y)))
p.lo <- c(1,2,5)/100
(P.p <- c(p.lo, (1:9)/10, p.hi))

## Proof that "quantile()" is just order statistics for this `Ns' :
iq.p <- 1 + P.p * (Ns - 1)
U <- runif(Ns)
stopifnot(identical(sort(U)[iq.p],  quantile(U, P.p, names=FALSE)))

P.dip <- matrix(nrow= length(nn), ncol = length(P.p),
		dimnames = list(as.character(nn), formatC(P.p,w=1,digits=7)))
dip.n <- numeric(Ns)
Tcpu <- 0
set.seed(963)
for(n in nn) {
    cat("n=",n,":")
    cpu <- system.time(for(i in 1:Ns)
                   {
                       if(i %% (Ns%/% 10) == 0) cat(i,"")
                       else if(i %% (Ns%/% 100) == 0) cat(".")
                       dip.n[i] <- dip(runif(n))
                   }
                       )[1:3]
    P.dip[paste(n),] <- quantile(dip.n, p = P.p, names=FALSE)
    cat("\nn=",n,", cpu=", paste(formatC(cpu), collapse=", "),"\n\n")
    Tcpu <- Tcpu + cpu
}

save(nn, Ns, P.p, P.dip, file="dipSim_1e5.rda", compress=TRUE)

cat("\nTotal CPU = ", paste(formatC(Tcpu), collapse=", "),"\n\n")

##--> File "dip-simul-sess"	to the CPU times .. !
q('no')
