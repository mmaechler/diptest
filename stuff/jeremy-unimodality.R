##-*- mode: R; kept-new-versions: 21; kept-old-versions: 12; -*-

#### From  http://www.stat.washington.edu/wxs/Stat593-s03/Code/jeremy-unimodality.R
####
#### After recommendation by Mark Difford,
####	on R-help Sep 09, 2008,  Subject  "Re: Modality Test"
#### https://stat.ethz.ch/pipermail/r-help/2008-September/173308.html

## TODO [MM] : compare this "seriously" with dip() from my package 'diptest'
##
####-----------------------------------------------------------------------------------
## Diagnostic plots for clustering and the DIP test for unimodality
## Code written by Jeremy Tantrum, Winter 2003
##=================================================================
##
## plot.ucdf(x) - plots the cdf of x and the closest unimodal cdf of x.
##
## plot.silverman(x) - plots the unimodal Gaussian smoother closest to the
##                     x and the closest bimodal Gaussian smoother.
##
## calcdip(x) - calculates the dip test statistic for x, using the mode found
##              by the closest unimodal Gaussian smoother.
##
## unisample(cd.out,n) - generates a sample from the unimodal distribution
##                       returned by the output of calcdip.
##
##
## An example of it working: Olive oil data - region 2 - are areas 5 and 6
## different.
## ---> see new file  ./jeremy-unimodality-olives.R
##                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##=================================================================

plot.ucdf <- function(x,plot.it = TRUE,COL1 = 1,COL2 = 2,LWD1 = 3,LWD2 = 2)
{
    x.cdf <- (1:length(x))/length(x)
    plot(c(min(x),sort(x)),c(0,x.cdf),type = 'n',xlab = "data",ylab = "cdf")
    lines(c(par('usr')[1],sort(c(x,x)),par('usr')[2]),sort(c(0,0,x.cdf,x.cdf)),
          lwd = LWD1,col = COL1)
    h.0 <- critwidth(x,4)
    x.f0 <- density(x,width = h.0$low)
    x.mode <- x.f0$x[order(x.f0$y)[length(x.f0$x)]]
    x.mode2 <- x[order(abs(x-x.mode))[1]]
    if(x.mode2 < sort(x)[4])
        x.mode2 <- sort(x)[4]
    if(x.mode2 > sort(x)[length(x)-4])
        x.mode2 <- sort(x)[length(x)-4]
    x.ord <- sort(x)

    x.split <- x.ord < x.mode2
    x.split2 <- x.ord >= x.mode2
    hull1 <- chull(c(x.ord[x.split],x.mode2),c(0,x.cdf[x.split]))
    n.1 <- sum(x.split)+1
    hlist <- c(hull1,hull1,hull1)
    start <- (1:length(hlist))[hlist == n.1][2]
    if(hlist[start+1] == 1) hlist <- rev(hlist)
    start <- (1:length(hlist))[hlist == n.1][2]
    x.hull1 <- n.1
    i <- start
    if(findslope(x.ord,x.cdf,1,hlist[i]) > findslope(x.ord,x.cdf,1,hlist[i+1]))
        while(hlist[i] > 1) {
            i <- i + 1
            x.hull1 <- c(x.hull1, hlist[i])
        }
    else
        while(hlist[i] > 1) {
            i <- i - 1
            x.hull1 <- c(x.hull1, hlist[i])
        }
    x.hull1 <- sort(x.hull1)
    hull2 <- chull(x.ord[x.split2],x.cdf[x.split2])
    n.2 <- sum(x.split2)
    hlist <- c(hull2,hull2,hull2)
    start <- (1:length(hlist))[hlist == n.2][2]
    if(hlist[start+1] == 1) hlist <- rev(hlist)
    start <- (1:length(hlist))[hlist == n.2][2]
    x.hull2 <- n.2
    i <- start
    if(findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i]) <
       findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i+1]))
        while(hlist[i] > 1) {
            i <- i +1
            x.hull2 <- c(x.hull2, hlist[i])
        }
    else
        while(hlist[i] > 1) {
            i <- i - 1
            x.hull2 <- c(x.hull2, hlist[i])
        }
    x.hull2 <- sort(x.hull2)
    lines(c(x.ord[x.split],x.mode2)[x.hull1],c(0,x.cdf[x.split])[x.hull1],
          col = COL2,lwd = LWD2)
    lines(x.ord[x.split2][x.hull2],x.cdf[x.split2][x.hull2],col = COL2,lwd = LWD2)

    hull.out <- list(hull1 = x.hull1, hull2 = x.hull2+n.1-1,
                     x = x.ord[sort(unique(c(x.hull1,x.hull2+n.1-1)))],
                     y = x.cdf[sort(unique(c(x.hull1,x.hull2+n.1-1)))])
    delta <- rep(0,length(x))
    for(i in 1:length(x))
        delta[i] <- abs(x.cdf[i] - fofx(sort(x)[i],hull.out))

    return(invisible(list(cdf = x.cdf,unicurve = hull.out, mode = x.mode2,
                          dip = max(delta), delta = delta,
                          dipwhere = order(delta)[length(delta)] )))
}

##-----------------------------------------------------------------

plot.silverman <- function(x,COL1 = 1,COL2 = "grey",...)
{
  h.0 <- critwidth(x,4,tol = 0.0001)$high
  h.1 <- critwidth2(x,h.0,tol = 0.001)$high
  f.c <- density(x,window = "g",width = h.0,n = 100)
  f.n <- density(x,width = h.1,n = 100)
  plot(c(f.c$x,f.n$x),c(f.c$y,f.n$y),type = 'n',xlab = "",ylab = "density",...)
  lines(f.n,lwd = 3,col = COL1)
  lines(f.c,lwd = 3,col = COL2)
  points(x,rep(0,length(x)),pch = '|')
}

## BELOW HERE ARE functions used by plot.silverman:

critwidth <- function(g.data,start,tol = 0.001,n.points = 200)
{
    if(is.unimodal(density(g.data,window = "g",width = start,n = n.points))) {
        high <- start
        low <- start/2
        while(is.unimodal(density(g.data,window = "g",width = low,n = n.points)))
            low <- low/2
    }
    else {
        low <- start
        high <- start*2
        while(!is.unimodal(density(g.data,window = "g",width = high,n = n.points)))
            high <- high*2
    }
    ## is.unimodal(low)=F and is.unimodal(high)=T
    while(high-low > tol) {
        wdth <- 0.5 * (high+low)
        if(is.unimodal(density(g.data,window = "g",width = wdth,n = n.points)))
            high <- wdth
        else
            low <- wdth
    }
    return(list(low = low,high = high))
}

###########################################################################
critwidth2 <- function(g.data,h.0,tol = 0.001,n.points = 200)
{
    ## h.0 is the critical width for a is.unimodal
    start <- h.0 + 2 * tol
    if(is.bimodal(density(g.data,window = "g",width = start,n = n.points))) {
        high <- start
        low <- start/2
        while(is.bimodal(density(g.data,window = "g",width = low,n = n.points)))
            low <- low/2
    }
    else {
        low <- start
        high <- start*2
        while(!is.bimodal(density(g.data,window = "g",width = high,n = n.points)))
            high <- high*2
    }
    ## is.unimodal(low)=F and is.unimodal(high)=T
    while(high-low > tol) {
        wdth <- 0.5 * (high+low)
        if(is.bimodal(density(g.data,window = "g",width = wdth,n = n.points)))
            high <- wdth
        else
            low <- wdth
    }
    return(list(low = low,high = high))
}

is.unimodal <- function(dens)
{
## dens is a list of dens$x and dens$y
  cdf <- cumsum(dens$y)
  n <- length(cdf)
  cdf.diff1 <- cdf[-1] - cdf[-n]
  cdf.diff2 <- cdf.diff1[-1] - cdf.diff1[-(n-1)]
  return(!any(order(-sign(cdf.diff2)) - 1:(n-2) > 0))
}

is.bimodal <- function(dens)
{
 ## dens is a list of dens$x and dens$y
  cdf <- cumsum(dens$y)
  n <- length(cdf)
  cdf.diff1 <- cdf[-1] - cdf[-n]
  cdf.diff2 <- cdf.diff1[-1] - cdf.diff1[-(n-1)]
  return(sum(sign(cdf.diff2)[-1] - sign(cdf.diff2)[-(n-2)] < 0) <= 2)
}

##-----------------------------------------------------------------

calcdip <- function(x, plot.it = TRUE, calc.it = TRUE)
{
    x <- sort(x)
    stopifnot((n <- length(x)) >= 4)
    h.0 <- critwidth(x,4)
    x.f0 <- density(x,width = h.0$low)
    x.mode <- x.f0$x[order(x.f0$y)[length(x.f0$x)]]
    x.mode2 <- x[order(abs(x-x.mode))[1]]
    if(x.mode2 < x[4])
        x.mode2 <- x[4]
    if(x.mode2 > x[n-4])
        x.mode2 <- x[n-4]
    x.cdf <- (1:n)/n
    hull.out <- findhulls(x,x.cdf,x.mode2,plot.it = plot.it,xlab = "",ylab = "CDF")
    delta <- rep(0,n)
    if(calc.it)
        for(i in 1:n)
            delta[i] <- abs(x.cdf[i] - fofx(x[i],hull.out))
    return(list(dip = max(delta),unicurve = hull.out))
}

## and the other functions needed:
unisample <- function(hull.out,size)
{
    n <- length(hull.out$x)
    min.x <- hull.out$x[1] - hull.out$y[1]/findslope(hull.out$x,hull.out$y,1,2)
    probs <- hull.out$y[-1] - c(0,hull.out$y[-c(1,n)])
    where <- sample(1:(n-1),size,replace = TRUE,prob = probs)
    out <- numeric(0)
    x <- c(min.x,hull.out$x[-1])
    for(i in 2:n) {
        x.s <- sum(where == i-1)
        if(x.s > 0)
            out <- c(out,runif(x.s,x[i-1],x[i]))
    }
    return(out)
}

findhulls <- function(x.ord,x.cdf,x.mode,plot.it = TRUE,...)
{
  x.split <- x.ord <= x.mode
  x.split2 <- x.ord >= x.mode
  hull1 <- chull(x.ord[x.split],x.cdf[x.split])
  n.1 <- sum(x.split)
  hlist <- c(hull1,hull1,hull1)
  start <- (1:length(hlist))[hlist == n.1][2]
  if(hlist[start+1] == 1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist == n.1][2]
  x.hull1 <- n.1
  i <- start
  if(findslope(x.ord,x.cdf,1,hlist[i]) > findslope(x.ord,x.cdf,1,hlist[i+1]))
    while(hlist[i] > 1) {
        i <- i + 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  else
    while(hlist[i] > 1) {
        i <- i - 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  x.hull1 <- sort(x.hull1)
  hull2 <- chull(x.ord[x.split2],x.cdf[x.split2])
  n.2 <- sum(x.split2)
  hlist <- c(hull2,hull2,hull2)
  start <- (1:length(hlist))[hlist == n.2][2]
  if(hlist[start+1] == 1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist == n.2][2]
  x.hull2 <- n.2
  i <- start
  if(findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i]) <
     findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i+1]))
    while(hlist[i] > 1) {
        i <- i +1
        x.hull2 <- c(x.hull2, hlist[i])
      }
  else
    while(hlist[i] > 1) {
        i <- i - 1
        x.hull2 <- c(x.hull2, hlist[i])
  }
  x.hull2 <- sort(x.hull2)
  if(plot.it) {
    plot(x.ord,x.cdf,...)
    lines(x.ord[x.split][x.hull1],x.cdf[x.split][x.hull1])
    lines(x.ord[x.split2][x.hull2],x.cdf[x.split2][x.hull2])
  }
  return(list(hull1 = x.hull1, hull2 = x.hull2+n.1 -1,
              x = x.ord[sort(unique(c(x.hull1,x.hull2+n.1-1)))],
              y = x.cdf[sort(unique(c(x.hull1,x.hull2+n.1-1)))]))
}

###########################################################################
findslope <- function(x,y,i,j)  return((y[j] - y[i])/(x[j]-x[i]))

###########################################################################
fofx <- function(x,hull.out)
{
    n <- length(hull.out$x)
    if(x <= hull.out$x[1])
        return(0)
    if(x >= hull.out$x[n])
        return(1)
    where <- (1:n)[order(c(x,hull.out$x)) == 1]
    return( (x-hull.out$x[where-1])/(hull.out$x[where]-hull.out$x[where-1]) *
           (hull.out$y[where]-hull.out$y[where-1]) + hull.out$y[where-1])
}

