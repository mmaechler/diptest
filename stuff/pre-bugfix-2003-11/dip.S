### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.S,v 1.3 1994/07/29 08:48:30 maechler Exp $

hartigan.data <- scan()
 30 33 35 36 37 37 39 39 39 39 39 40 40 40 40 41 42 43 43 43 44 44 45 45 46
 46 47 47 48 48 48 49 50 50 51 52 52 53 53 53 53 53 54 54 57 57 59 60 60 60
 61 61 61 61 62 62 62 62 63 66 70 72 72

dip <- function(x, full.result = FALSE)
{
  ## Purpose: Compute the "Dip test for unimodality" (statistic)
  ## -------------------------------------------------------------------------
  ## Arguments: x: the Data, full.result: return also 'xl', 'xu' (modal interv.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler <maechler@stat.math.ethz.ch>, Jul 94
  ##      from 1st version of   Dario Ringach <dario@wotan.cns.nyu.edu>
  
  if(!is.loaded(symbol.C("diptst"))) dyn.load2("dip.o")

  n <- length(x)
  ii <- if(full.result) 1:10 else  "dip"
  simplify <- function(lis) if(length(lis)==1) lis[[1]] else lis
  simplify(.C("diptst",
	      x   = as.single(sort(x)),
	      n   = as.integer(n),
	      dip = single(1),
	      xl  = single(1),
	      xu  = single(1),
	      ifault= integer(1),
	      gcm =   integer(n),
	      lcm =   integer(n),
	      mn  =   integer(n),
	      mj  =   integer(n))[ii])
  ##-   if(z$ifault!=0)  #-- something not ok, but this is IMPOSSIBLE here
  ##-     stop(paste("Problem -- C 'message' : ifault = ", z$ifault))
}

dip(hstart)## 0.03998769
dip(ship)  ## 0.02328392
dip(iris)  ## 0.02874999
dip(hartigan.data) ## Should give   0.05952381
str(dip.hh <- dip(hartigan.data, full = T))

## NOTA BENE:  dip(1:10) gives INFINITE LOOP (error in Hartigan's Algorithm!)

dip.1.10_rep(0,100)
unix.time(for(i in 1:100) dip.1.10[i] <- dip(jitter(1:10)))
##--> 2.3 CPU sec. on Sparstation 10 (ingrid)
summary(dip.1.10)
##-  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##-  0.05 0.05757 0.06016 0.05979 0.06237 0.06686
##-  0.05   0.058 0.06027 0.06011 0.06231 0.06696
##-- ~~~~ EXAKT (9 times out of 200)

cdf <- function(x, X=x)
{
  ## Purpose: Empirical distribution function F_{X}(x)
  ## -------------------------------------------------------------------------
  ## Arguments: x: vector of ARGuments;  X : data
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 28 Jul 94, 18:03
  apply(outer(X,x,"<="),2,sum) / length(X)
}
## For plotting,  plot.step(X)  is much better !



dip0 <- function(x)
{
  ## Load in ANY case ! if(!is.loaded(symbol.C("diptst")))
  dyn.load2("dip0.o", userlibs="-lf2c")
  n <- length(x)
  .C("diptst",
     as.single(sort(x)),
     as.integer(n),
     ans = single(1),
     as.single(0),
     as.single(0),
     as.integer(0),
    integer(n),
    integer(n),
    integer(n),
    integer(n)) $ans
}
dip0(hartigan.data) ##-> [1] 0.05952381

dipF <- function(x)
{ 
  if(!is.loaded(symbol.For("diptst")))
  dyn.load2("dipF.o")
  n <- length(x)
  .Fortran("diptst",
     as.single(sort(x)),
     as.integer(n),
     ans = single(1),
     as.single(0),
     as.single(0),
     as.integer(0),
    integer(n),
    integer(n),
    integer(n),
    integer(n)) $ans
}
dipF(hartigan.data) ##-> [1] 0.05952381

### Test if  dip.c  is still okay, after all my hacking at the code
###                 (dipF.f  which is used by dipF(..) is ORIGINAL !)

for(i in 1:100){x_runif(50); if(dip(x) != dipF(x)) { xx<<-x; stop("DIFFERENT")}}
for(i in 1:100) { 
  x_c(rnorm(100), 5+ rnorm(20))
  if(dip(x) != dipF(x)) { xx<<-x; stop("DIFFERENT") } }
