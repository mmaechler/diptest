### From the large (Ns = 100001) simulation in
### ./new-simul1e5.R
###  ~~~~~~~~~~~~~~~
# "wrong" dip stat: load("/u/maechler/R/Pkgs/diptest/stuff/dipSim_1e5.rda")
load("/u/maechler/R/Pkgs/diptest/stuff/dipSim_1e6.rda")

stopifnot(identical(P.p, as.numeric(colnames(P.dip))),
          identical(nn,  as.numeric(rownames(P.dip))))

names(dimnames(P.dip)) <- c("n","Pr")
## new data set!
qDiptab <- P.dip

      Pp <- P.p  [  P.p < 1]# not max() = 100% percentile
P.dip.Rn <- P.dip[, P.p < 1]*sqrt(nn)
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

matplot(nn, P.dip.Rn, type='n', xlab = 'n  [log scale]', ylab = y.tit,
        xlim = range(1.5,nn), log='x', main= nqdip.tit)
mPerclegend(1.0, 0.95, Pp)
matlines(nn, P.dip.Rn[, "0.5"], lwd = 3, col = "dark gray")
matlines(nn, P.dip.Rn, type = 'o')
mtext(paste(Ns, " simulated samples"), 3, line=0)
mtext(paste("/u/maechler/R/Pkgs/diptest/stuff/", date(),sep="\n"),
      4, cex=.8, adj=0)

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
matplot(nn, P.dip.Rn[, smP], type='n', log='xy', xlim = range(1.5,nn),
        xlab = 'n  [log scale]',ylab = y.titL, main = nqdip.tit)
mPerclegend(1.0, .75, Pp[smP])
matlines(nn, P.dip.Rn[, "0.5"], lwd = 3, col = "dark gray")
matlines(nn, P.dip.Rn[, smP], type = 'o')
subtit(paste(Ns, " simulated samples"))



### only larger N   to see if it became constant:
nP <- sum(nL <- nn > 100)
matplot(nn[nL], P.dip.Rn[nL,], type='o', xlim = range(50,nn[nL]),
        log='x', xlab = 'n  [log scale]', ylab = y.tit, main = nqdip.tit)
matlines(nn[nL], P.dip.Rn[nL, "0.5"], lwd = 3, col = "dark gray")
##matlines(nn[nL], P.dip.Rn[nL,], type = 'o')
## FIXME: col, lty are WRONG
legend(45, 0.95, legend=
       rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
                 paste(100*Pp,"%",sep=""), sep=": ")),
       col = rev(rep(1:6,length=nP)), lty = rev(rep(1:5,length=nP)),bty='n')


mult.fig(9, main = "dip(U[0,1]) distribution {simulated} -- for small n")
for (cn in nn[1:9]) {
    plot(P.dip[paste(cn),], P.p,
         xlab = "dip = d(x[1 .. n])", ylab = expression(P(D >= d)),
         type = 'o', cex = 0.6, main = paste("n = ",cn))
    abline(h=0:1, col="gray")
}
###

##-- now really
pDip <- function(D, n)
{
  ## Purpose:  compute the P-value of dip(X) = D  under H0: X ~ U[0,1]
  ## ----------------------------------------------------------------------
  ## Arguments: D : result of dip();  n : sample size
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 Jul 2003, 09:39

}

## Another thing: Draw "density" from quantiles only:
## Use derivative of cubic spline interpolation or so?


P2dens <- function(x, probs, eps.p = 1e-7, xlim = NULL, f.lim = 0.5)
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
        stop("`x' and `probs' must have same length")
    if(any(0 > probs | probs > 1))
        stop("`probs' must be in [0,1]")
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

    library(splines) ## even better: use constrained splines: F monotone
    Fspl <- interpSpline(x, probs)
    rm(x,probs)
    function(x) predict(Fspl, x, deriv = 1) $y
}

d1 <- P2dens(-3:3, pr=pnorm(-3:3))
plot(d1, -7, 7, n = 501)
# quite fine :
curve(dnorm, col = 2, add=TRUE, n = 501)

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
dDips <- apply(P.dip, 1, function(qd) P2dens(qd, pr = P.p))
## Error in interpSpline..(.) : values of x must be distinct << !
