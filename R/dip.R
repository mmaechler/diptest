### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.S,v 1.4 1994/07/29 10:01:12 maechler Exp $

dip <- function(x, full.result = FALSE, debug = FALSE)
{
  n <- length(x)
  (function(lis) if(length(lis)==1) lis[[1]] else lis
   ) (.C("diptst",
         x   = as.double(sort(x)),
         n   = as.integer(n),
         dip = double(1),
         xl  = double(1),
         xu  = double(1),
         ifault= integer(1),
         gcm =   integer(n),
         lcm =   integer(n),
         mn  =   integer(n),
         mj  =   integer(n),
         debug= as.logical(debug))[if(full.result) TRUE else "dip"])
##-   if(z$ifault)  #-- something not ok, but this is IMPOSSIBLE here
##-     stop(paste("Problem -- C 'message' : ifault = ", z$ifault))
}
