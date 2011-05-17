### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

###-- $Id: dip.R,v 1.8 2011/05/16 06:34:19 maechler Exp maechler $

dip <- function(x, full.result = FALSE, min.is.0 = FALSE, debug = FALSE)
{
  if(full.result) cl <- match.call()

  n <- as.integer(length(x))
  x <- sort(unname(x), method="quick")
  ## be careful to "duplicate" (as have DUP=FALSE):
  min.is.0 <- as.logical(min.is.0)
  debug <- as.integer(debug)# FALSE/TRUE or 2, 3, ...
  r <- .C("diptst",
          x   = as.double(x),
          n   = n,
          dip = double(1),
          lo.hi = integer(4),
          ifault= integer(1),
          gcm =   integer(n),
          lcm =   integer(n),
          mn  =   integer(n),
          mj  =   integer(n),
          min.is.0 = min.is.0,
          debug = debug,
          DUP = FALSE,
          PACKAGE = "diptest")[if(full.result) TRUE else "dip"]
  ## if(r$ifault)	#-- something not ok, but this is IMPOSSIBLE here
  ##     stop("Problem -- should not happen, please report! --  ifault = ",
  ##          r$ifault)
  if(full.result) {
      l.GL <- r$lo.hi[3:4]
      length(r$gcm) <- l.GL[1]
      length(r$lcm) <- l.GL[2]
      length(r$lo.hi) <- 2L
      u <- x[r$lo.hi]
      structure(class = "dip",
		c(list(call = cl), r, list(xl = u[1], xu = u[2])))
  }
  else r[[1]]
}

print.dip <- function(x, digits = getOption("digits"), ...)
{
    stopifnot(is.integer(n <- x$n), is.numeric(D <- x$dip),
              length(lh <- x$lo.hi) == 2)
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("n = ", n,".  Dip statistic, D_n = ",
        format(D, digits=digits),"; 2n * D_n = ",
        format(2*n* D, digits=digits),
        sprintf("\n Modal interval [x[%d], x[%d]] = [%g, %g]\n",
                lh[1], lh[2], x$x[lh[1]], x$x[lh[2]]),
        sprintf(" GCM and LCM have %d and %d nodes, respectively\n",
                sum(x$gcm != 0), sum(x$lcm != 0)),
        sep="")
    invisible(x)
}

aLine <- function(r.dip, lType = c("gcm","lcm"),
                  type = "b", col="red3", lwd=1.5)
{
    lType <- match.arg(lType)
    stopifnot(is.numeric(x <- r.dip$x),
              is.integer(n <- r.dip$n),
              is.integer(i <- r.dip[[lType]]) # 'gcm' or 'lcm' or component
              )
    e <- if(lType =="gcm") .01*min(diff(unique(x))) else 0
    i <- i[i != 0]
    lines(x[i], ecdf(x)(x[i] - e),
          type=type, col=col, lwd=lwd)
}

plot.dip <- function(x, colG="red3", colL="blue3", colM="forest green",
                     ## for plot.stepfun():
                     do.points=(n < 1000), col.points=par("col"), col.hor=col.points,
                     doModal=TRUE, doLegend=TRUE, ...)
{
    stopifnot(is.integer(n <- x$n), is.numeric(D <- x$dip),
              length(lh <- x$lo.hi) == 2)
    Fn <- ecdf(x$x)
    ## and now manipulate the call such that it's plotted nicely
    cl <- x$call[1:2]
    cl[[1]] <- as.name("ecdf") ; names(cl)[2] <- ""
    attr(Fn, "call") <- cl
    chD <- formatC(D, digits=pmax(3, getOption("digits")-2))
    tit <- bquote("Dip" ~~ {D[n] == D[.(n)]} == .(chD))
    plot(Fn, do.points=do.points, col.points=col.points, col.hor=col.hor,
         verticals=TRUE, col.vert = "sky blue", lwd=2, ...)
    title(tit, adj = 0, line = 1.25)
    aLine(x, "gcm", col=colG)
    aLine(x, "lcm", col=colL)
    if(doModal) {
        x12 <- x$x[lh]
        abline(v= x12, col = colM, lty = 2)
        axis(3, at=x12, label = expression(x[L], x[U]),
             tick=FALSE, padj = .9, col.axis = colM)
    }
    if(doLegend)
        legend("topleft", bty = "n",
               c("greatest convex minorant GCM",
                 "least concave majorant    LCM"),
               col =c(colG, colL), lty=1, lwd=1.5)
    invisible()
}

