### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.S,v 1.4 1994/07/29 10:01:12 maechler Exp $

###--- 1st page of this file should be  source(.)able !! -----------

if(!interactive()) {
  options(echo = T) #-- when being run as batch [FAILS for source() !!]
  print(date())
}


DIP.dir <- getenv("PWD")

hartigan.data <-
  c(30,33,35,36,37,37,39,39,39,39,39,40,40,40,40,41,42,43,43,43,44,44,45,45,46,
    46,47,47,48,48,48,49,50,50,51,52,52,53,53,53,53,53,54,54,57,57,59,60,60,60,
    61,61,61,61,62,62,62,62,63,66,70,72,72)

dip <- function(x, full.result = FALSE, debug = F)
{
  ## Purpose: Compute the "Dip test for unimodality" (statistic)
  ## -------------------------------------------------------------------------
  ## Arguments: x: the Data, full.result: return also 'xl', 'xu' (modal interv.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler <maechler@stat.math.ethz.ch>, Jul 94
  ##      from 1st version of   Dario Ringach <dario@wotan.cns.nyu.edu>
  
  if(!is.loaded(symbol.C("diptst"))) 
    dyn.load(paste(DIP.dir,"dip.o", sep="/"))

  n <- length(x)

  ii <- if(full.result) 1:11 else  "dip"
  simplify <- function(lis) if(length(lis)==1) as.single(lis[[1]]) else lis
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
	      mj  =   integer(n),
	     debug= as.integer(debug))[ii])
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

p.cdf <- function(X, type='b',...)
{
  ## Purpose:  PLOT empirical distribution function
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 30 Jul 94, 11:41
  X _ sort(X)
  plot(X, cdf(X), type=type, ...)
}

## NOTE: For plotting, Martin Maechler's   plot.step(X)  is MUCH better !


dip0 <- function(x, debug=F)
{
  ## Load in ANY case ! if(!is.loaded(symbol.C("diptst")))
  dyn.load(paste(DIP.dir,"dip0.o", sep="/"))
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
    integer(n),
     as.integer(debug)) $ans
}
dip0(hartigan.data) ##-> [1] 0.05952381

## Load my current "diptst" C-symbol !
dyn.load(paste(DIP.dir,"dip.o", sep="/"))


dipF <- function(x)
{ 
  if(!is.loaded(symbol.For("diptst")))
    dyn.load(paste(DIP.dir,"dipF.o", sep="/"))
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


##-- NOTE: for small 'n' (n=4,5, even 8, 10 (rarely!)) sometimes INFINITE loop!!

###===> FIXED in big.c  !!

x4.ok <- x4.bad <- NULL
xx4_round(runif(4),3); cat(format(dipF(xx4))); x4.ok_cbind(x4.ok,xx4)
##-- C-c (kill) the above line if running long;  then do
x4.bad <- cbind(x4.bad, xx4, deparse.level=0)

##-- After a while (MANUALLY !)
dim(x4.bad) #  4 10
dim(x4.ok)  #  4 34


##-- n=10 even can be bad: 1st example:  1:10
## 2nd example: 
n_10; set.seed(99);for(i in 1:3635) x _ runif(n)
x10.bad _ round(sort(x),4)
dipF(round(x10.bad,3)) #-- .05 (minimal value)
plot.step(x10.bad);  subtit(vcat(x10.bad,sep=",  "))

u.dev.default() 
nb_ncol(x4.bad); mult.fig(nb); for(i in 1:nb) plot.step(x4.bad[,i],cad.lag=F,main="")
u.dev.default() 
no_ncol(x4.ok);  mult.fig(no); for(i in 1:no) plot.step(x4.ok[,i],cad.lag=F,main="")
