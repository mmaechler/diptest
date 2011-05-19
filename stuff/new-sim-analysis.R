### From the large (Ns = 100001) simulation in
### ./new-simul.R
###  ~~~~~~~~~~~~

Ns <- 1000001

if(require("diptest") && is.character(data(qDiptab))) {
    dn <- dimnames(qDiptab)

    nn  <- as.numeric(dn[[1]])
    P.p <- as.numeric(dn[[2]])

} else { ## no longer needed, now we have  'qDiptab'
    setwd("/u/maechler/R/Pkgs/diptest/stuff")

    ## "wrong" dip stat: load("/u/maechler/R/Pkgs/diptest/stuff/dipSim_1e5.rda")
    load("dipSim_1e6.rda")

    stopifnot(identical(P.p, as.numeric(colnames(P.dip))),
              identical(nn,  as.numeric(rownames(P.dip))))

    names(dimnames(P.dip)) <- c("n","Pr")

    data(qDiptab, package="diptest")
    identical(P.dip, qDiptab) # !

    attr(P.dip, "N_1") <- as.integer(Ns - 1)

    ## new data set!
    qDiptab <- P.dip

} # end {else: no longer needed}


Pp <- P.p  [  P.p < 1]# not max() = 100% percentile
qDip.Rn <- qDiptab[, P.p < 1]*sqrt(nn)
nP <- length(Pp)

## Titles
nqdip.tit <-
    expression(sqrt(n) %*% " percentile of  "* list(dip(X),
        " "* X *" ~ " * U * group("[",list(0,1),"]"))*
        "   vs.  " * n)
y.tit <- expression(sqrt(n) *" " %*% " " * qDip)
(y.titL <- substitute(F * "  [log scaled]", list(F = y.tit[[1]])))

mPerclegend <- function(x,y, pr) {
    nP <- length(pr)
    legend(x, y, legend=
           rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
                     paste(100*pr,"%",sep=""), sep=": ")),
           col = rev(rep(1:6,length=nP)),
           lty = rev(rep(1:5,length=nP)), bty='n')
}

if(FALSE)## for paper --- write a vignette !! ---
sfsmisc::pdf.latex("dip_critical.pdf")
## for print out, or "sending along":

pdf.do("dip_critical.pdf", paper = "a4")

matplot(nn, qDip.Rn, type = "n", yaxt = "n",
        xlab = 'n  [log scale]', ylab = y.tit,
        xlim = range(1.5,nn), log='x', main= nqdip.tit)
abline(h = seq(0.1, 0.9, by = 0.05), col = "gray80", lty = 3)
op <- par(las=2); for(i in c(2,4)) axis(i, at = seq(0.1,0.9, by=0.1)); par(op)
mPerclegend(1.0, 0.95, Pp)
matlines(nn, qDip.Rn[, "0.5"], lwd = 3, col = "dark gray")
matlines(nn, qDip.Rn, type = 'o')
mtext(paste(Ns, " simulated samples"), 3, line=0)
mtext("© Martin Maechler, ETH Zurich", 1, line = 3.2, adj = 1)

pdf.end()

mtext(paste(getwd(), date(), sep="\n"), 4, cex=.8, adj=0)

## "research" :
## 1. prove that min(dip) = 1/(2 * n)
##	{I'm sure this follows from the "string" equivalence;
##       but that is not according to theory D(F) dip definition, where
##	 D(F) = 0  <==> F itself is unimodal
##    Can we prove that  dip(F_n) = 1/(2n)  ===>  F_n is unimodal ?
##
## 2. derive the function  ff(n) :=  Prob[dip = 1/(2 * n)]


## log y: --->> more symmetric distributions, but still skewed to the right
## -----        (asymptotic) is of interest, but also: how to do interpolation
nP <- sum(smP <- 0.01 < Pp & Pp <= .999)# only "relevant subset"
matplot(nn, qDip.Rn[, smP], type='n', log='xy', xlim = range(1.5,nn),
        xlab = 'n  [log scale]',ylab = y.titL, main = nqdip.tit)
mPerclegend(1.0, .75, Pp[smP])
matlines(nn, qDip.Rn[, "0.5"], lwd = 3, col = "dark gray")
matlines(nn, qDip.Rn[, smP], type = 'o')
mtext(paste(Ns, " simulated samples"), side = 3, line = 0)



### only larger N   to see if it became constant:
nN <- sum(nL <- nn > 100)
matplot(nn[nL], qDip.Rn[nL,], type='o', xlim = range(50,nn[nL]),
        log='x', xlab = 'n  [log scale]', ylab = y.tit, main = nqdip.tit)
mtext(paste(Ns, " simulated samples"), side = 3, line = 0)
mPerclegend("topleft", NULL, Pp)
matlines(nn[nL], qDip.Rn[nL, "0.5"], lwd = 4,
         col = adjustcolor("black", 0.4))


##-- Asymptotic : see more in ./asymp-distrib.R
##                              ~~~~~~~~~~~~~~~


mult.fig(9, main = "dip(U[0,1]) distribution {simulated} -- for small n")
for (cn in nn[1:9]) {
    plot(qDiptab[paste(cn),], P.p,
         xlab = "dip = d(x[1 .. n])", ylab = expression(P(D >= d)),
         type = 'o', cex = 0.6, main = paste("n = ",cn))
    abline(h=0:1, col="gray")
}
###


## Another thing: Draw "density" from quantiles only:
## Use derivative of cubic spline interpolation or so?


P2dens <- function(x, probs, eps.p = 1e-7, xlim = NULL, f.lim = 0.5,
                   method = c("interpSpline", "monoH.FC", "natural"))
{
    ## Purpose: Density from probabilty/quantiles -- return a *function*
    ## ----------------------------------------------------------------------
    ## Arguments:     x: quantiles , x[1:n]
    ##            probs: probabilities, i.e., Pr{X <= x[i]} == probs[i]
    ##            eps.p: small value used in:
    ##            xlim: if(probs[1] > eps.p and/or probs[n] < 1-eps.p),  use
    ##                   x[0] = xlim[1] and/or x[n+1] = xlim[2]  with
    ##                   probs[0] = 0 and/or probs[n+1] = 1
    ##                  Per default, xlim = range(x) +  f.lim * c(-d,d),
    ##                  where d = diff(range(x)) = max(x) - min(x)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 14 Jul 2003, 14:57
    if((n <- length(x)) != length(probs))
        stop("'x' and 'probs' must have same length")
    if(any(0 > probs | probs > 1))
        stop("'probs' must be in [0,1]")
    if(is.unsorted(probs)) {
        sp <- sort(probs, index.return=TRUE)
        x <- x[sp$ix]
        probs <- sp$x
        rm(sp)
    }

    if((L <- probs[1] > eps.p) |
       (R <- probs[n] < 1-eps.p)) {
        d <- diff(r <- range(x))
        if(L) { x <- c(r[1] - f.lim*d, x); probs <- c(0, probs) }
        if(R) { x <- c(x, r[2] + f.lim*d); probs <- c(probs, 1) }
        rm(r,d,L,R)
    }

    method <- match.arg(method)
    switch(method,
	   "interpSpline" = {
	       library(splines)
	       Fspl <- interpSpline(x, probs)
	       rm(x,probs)
	       function(x) predict(Fspl, x, deriv = 1) $y
	   },
	   "monoH.FC" = {
	       f <- splinefun(x, probs, method="mono")
	       formals(f)[["deriv"]] <- 1
	       f
	   },
	   "natural" = {
	       f <- splinefun(x, probs, method="natural")
	       formals(f)[["deriv"]] <- 1
	       f
	   },
	   ## otherwise
	   stop("invalid method ", method))
}

d1 <- P2dens(-3:3, pr=pnorm(-3:3))
plot(d1, -7, 7, n = 501)
points(-3:3, d1(-3:3))
## quite fine :
curve(dnorm, col = 2, add=TRUE, n = 501)

## Now with *monotone* Hermite interpolation --- this is *WORSE* !!
d2 <- P2dens(-3:3, pr=pnorm(-3:3), method = "mono")
plot(d2, -7, 7, n = 501, add=TRUE, col="midnightblue")
points(-3:3, d2(-3:3))

## and more experiments suggest the best (here!) solution being "natural"
d2 <- P2dens(-3:3, pr=pnorm(-3:3), method = "natural")
plot(d2, -7, 7, n = 501, add=TRUE, col="midnightblue")
points(-3:3, d2(-3:3))
## well,
x <- seq(-7,7, len=1001)
all.equal(d1(x), d2(x), tol = 1e-15)# TRUE

str(get("Fspl", envir= environment(d1)))
##- List of 2
##-  $ knots       : num [1:9] -6 -3 -2 -1 0 1 2 3 6
##-  $ coefficients: num [1:9, 1:4] 0.00000 0.00135 0.02275 0.15866 0.50000 ...
##-  - attr(*, ...........

x0 <- c(0, 2^(-3:6))
d2 <- P2dens(x0, pr=pgamma(x0, shape = 1.5))
plot(d2, 0, max(x0), n = 501)
rug(x0)
# not bad, too :
curve(dgamma(x,shape=1.5), col = 2, add=TRUE, n = 501)


##-- But the thing I wanted fails because of point masses (left border):
dDips <- apply(qDiptab, 1, function(qd) P2dens(qd, pr = P.p))
## Error in interpSpline..(.) : values of x must be distinct << !
