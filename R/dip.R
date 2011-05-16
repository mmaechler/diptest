### S-interface to Hartigan's algorithm for "The dip test for unimodality"-- $Id$
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.S,v 1.4 1994/07/29 10:01:12 maechler Exp $

dip <- function(x, full.result = FALSE, dip.min.0 = FALSE, debug = FALSE)
{
  n <- length(x)
  x <- sort(unname(x), method="quick")
  ## be careful to "duplicate" (as have DUP=FALSE):
  min.is.0 <- as.logical(dip.min.0)
  debug <- as.logical(debug)
  r <- .C("diptst",
          x   = as.double(x),
          n   = n,
          dip = double(1),
          lo.hi = integer(2),
          ifault= integer(1),
          gcm =   integer(n),
          lcm =   integer(n),
          mn  =   integer(n),
          mj  =   integer(n),
          min.is.0 = min.is.0,
          debug = debug,
          DUP = FALSE,
          PACKAGE = "diptest")[if(full.result) TRUE else "dip"]
  ##-   if(z$ifault)  #-- something not ok, but this is IMPOSSIBLE here
  ##-     stop(paste("Problem -- C 'message' : ifault = ", z$ifault))
  if(full.result) c(r, {u <- x[r$lo.hi]; list(xl = u[1], xu = u[2])})
  else r[[1]]
}
