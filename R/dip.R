### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.S,v 1.4 1994/07/29 10:01:12 maechler Exp $

dip <- function(x, full.result = FALSE, debug = FALSE)
{
  ## Purpose: Compute the "Dip test for unimodality" (statistic)
  ## -------------------------------------------------------------------------
  ## Arguments: x: the Data, full.result: return also 'xl', 'xu' (modal interv.)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler <maechler@stat.math.ethz.ch>, Jul 94
  ##      from 1st version of   Dario Ringach <dario@wotan.cns.nyu.edu>

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
