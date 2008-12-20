###----- More  "full simulations" of the asymptotic limit:

setwd("/u/maechler/R/Pkgs/diptest/stuff")
## These all have  n.sim = 1000'001  samples of  dip(runif(N)):
load("dip8k.rda")
load("dip12k.rda")
load("dip20k.rda")

## And the empirical densities of the sqrt(n)-blown up look
## "practically" identical:
plot (density(sqrt(20000)* dip20k))
lines(density(sqrt(12000)* dip12k), col=2)
lines(density(sqrt(8000) * dip8k), col=3)

## Could Log-normal fit ?
plot (density(log(sqrt(20000)* dip20k)))
lines(density(log(sqrt(12000)* dip12k)), col=2)
lines(density(log(sqrt(8000) * dip8k)), col=3)
## looks quite symmetric (but is not quite, see below!)

x <- sort(sqrt(20000)* dip20k) # sorting, just for possible convenience
lx <- log(x)

if(FALSE) ## Attention! Plotting 1 mio points on Xlib("cairo") is deadly!
qqnorm(lx) #-> already clear that (log)normal does  NOT fit

library(MASS)
fd <- fitdistr(x, "lognormal")
fd
##       meanlog          sdlog
##   -0.9656529434    0.2074114508
##  ( 0.0002074113) ( 0.0001466620)
logLik(fd) # 1119766

fd. <- fitdistr(lx, "normal")
fd. # exact same parameters as above
logLik(fd.) # 154112.4 (df=2) --- very different to  fd's

dlnormFit <- function(x) do.call(dlnorm, c(list(x=x), coef(fd)))
dnormFit <- function(x) do.call(dnorm, c(list(x=x), coef(fd.)))
curve(dnormFit, add = TRUE, col = "tomato")
##-->  log-normal clearly does *not* fit !


fdg <- fitdistr(x, "gamma") ## this takes a few seconds!
fdg
##       shape         rate
##   23.09376812   59.34742081
##  ( 0.03242641) ( 0.08424100)
dgammaFit <- function(x) do.call(dgamma, c(list(x=x), coef(fdg)))

fdw <- fitdistr(x, "weibull") ## this takes a few seconds!
fdw <- fitdistr(x, "weibull", control=list(trace=2))
fdw
##       shape          scale
##   4.665147e+00   4.231184e-01
##  (3.277802e-03) (9.629305e-05)
dweibullFit <- function(x) do.call(dweibull, c(list(x=x), coef(fdw)))

## In original scale :
plot (density(sqrt(20000)* dip20k))
lines(density(sqrt(12000)* dip12k), col=2)
lines(density(sqrt(8000) * dip8k), col=3)
curve(dlnormFit, add = TRUE, col = "tomato")
## does not fit
curve(dgammaFit, add = TRUE, col = "blue3")
## is even worse than log-normal
curve(dweibullFit, add = TRUE, col = "forest green")
## is much worse even


## How does the upper tail look?
## log(1 - P) = log(P{X >= x})
