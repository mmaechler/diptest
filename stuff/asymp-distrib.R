###----- More  "full simulations" of the asymptotic limit:

setwd("/u/maechler/R/Pkgs/diptest/stuff")

## These all have  n.sim = 1000'001  samples of
## dip(runif(N)) for "large N"
## where produced by scripts such as  ./d20k_do.R
N.k.set <- c(8,12,16,20,24,32,36,40)
## or automatically
patt <- "^dip(.*)k\\.rda"
N.k.set <- sort(as.integer(sub(patt, "\\1", list.files(pattern = patt))))
## for now:
if(5 %in% N.k.set) { N.k.set <- N.k.set[N.k.set != 5] }

dip.nm <- function(N.k, file=FALSE)
    paste("dip", N.k, if(file)"k.rda" else "k", sep='')
for(N.k in N.k.set)
    load(dip.nm(print(N.k), file=TRUE))


d.dip <- function(N.k, scaleUp = TRUE)
{
    ## Purpose: Simulation data for N = N.k * 1000, possibly sqrt(N) scaled
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15 Apr 2009, 10:36
    N <- N.k * 1000
    nm <- dip.nm(N.k)
    if(scaleUp) sqrt(N) * get(nm) else get(nm)
}
Nk2char <- function(Nk) {
    dput(Nk, textConnection(".val", "w", local=TRUE))
    .val
}
mFormat <- function(n) sub("([0-9]{3})$", "'\\1", format(n))

names(N.k.set) <- paste(format(N.k.set), "'000", sep='')
t(sums <- sapply(N.k.set, function(Nk) summary(d.dip(Nk))))
##          Min. 1st Qu. Median   Mean 3rd Qu.   Max.
##  6'000 0.1777  0.3274 0.3759 0.3876  0.4357 0.9858
##  8'000 0.1795  0.3278 0.3765 0.3881  0.4364 0.9772
## 10'000 0.1750  0.3280 0.3767 0.3884  0.4368 0.9522
## 12'000 0.1779  0.3283 0.3770 0.3887  0.4370 0.9421
## 16'000 0.1768  0.3285 0.3774 0.3890  0.4373 1.0080
## 20'000 0.1689  0.3288 0.3774 0.3891  0.4373 0.9699
## 24'000 0.1814  0.3288 0.3775 0.3893  0.4375 1.0570
## 32'000 0.1758  0.3290 0.3779 0.3896  0.4380 0.9757
## 36'000 0.1724  0.3295 0.3783 0.3899  0.4382 1.0500
## 40'000 0.1705  0.3294 0.3781 0.3899  0.4382 1.0270
## 44'000 0.1742  0.3295 0.3782 0.3898  0.4381 0.9693
## 52'000 0.1752  0.3295 0.3783 0.3900  0.4383 0.9463
## 60'000 0.1755  0.3295 0.3785 0.3901  0.4384 1.0060
## 72'000 0.1729  0.3296 0.3785 0.3902  0.4386 1.0100

## And the empirical densities of the sqrt(n)-blown up look
## "practically" identical:
plot (density(d.dip(36)))
for(N.k in head(N.k.set,-1))
    lines(density(d.dip(N.k)))
mtext(paste("superimposed dip-distrib. for  N = 1000 * ", Nk2char(N.k.set)))

## Could Log-normal fit ?
plot (d <- density(log(d.dip(max(N.k.set)))))
for(N.k in head(N.k.set,-1))
    lines(density(log(d.dip(N.k)), bw = d$bw))
mtext(paste("superimposed log(dip)-distrib. for  N = 1000 * ",
            Nk2char(N.k.set)))

## Now only for larger N:
plot (d <- density(log(d.dip(max(N.k.set)))))
for(N.k in tail(head(N.k.set,-1), 3))
    lines(density(log(d.dip(N.k)), bw = d$bw))
mtext(paste("superimposed log(dip)-distrib. for N = 1000 * ",
            Nk2char(tail(N.k.set, 4))))


## looks quite symmetric (but is not quite, see below!)

x <- sort(sqrt(20000)* dip20k) # sorting, just for possible convenience
lx <- log(x)

if(FALSE) { ## Attention! Plotting 1 mio points on Xlib("cairo") is deadly!
    x11(type = "Xlib")
    qqnorm(lx) #-> already clear that (log)normal does  NOT fit
}
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
mtext("dnorm(<fit>) super-imposed", col = "tomato")
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
plot (density(d.dip(36)))
lines(density(d.dip(32)), col="gray")
lines(density(d.dip(24)), col="gray30")

curve(dlnormFit, add = TRUE, col = "tomato")
## does not fit
curve(dgammaFit, add = TRUE, col = "blue3")
## is even worse than log-normal
curve(dweibullFit, add = TRUE, col = "forest green")
## is much worse even

###---------> "back" to look at  CDFs ---------------------------------

Fn12k <- ecdf(d.dip(12))
Fn20k <- ecdf(d.dip(20))
Fn24k <- ecdf(d.dip(24))
Fn32k <- ecdf(d.dip(32))
Fn36k <- ecdf(d.dip(36))


plot(Fn20k, do.points = FALSE)# still slow
lines(ecdf(sqrt(12000) * dip12k), col=2, do.points = FALSE)

ks.test(sqrt(12000) * dip12k,
        sqrt(20000) * dip20k)
## 	Two-sample Kolmogorov-Smirnov test

## data:  sqrt(12000) * dip12k and sqrt(20000) * dip20k
## D = 0.0029, p-value = 0.0004664  <<<--- ***
## alternative hypothesis: two-sided

## and comparing the 8000 with 12'000 case even
## has  p-value = 4.98 e-5

## Now, with new "24'000":
ks.test(sqrt(24000) * dip24k,
        sqrt(20000) * dip20k)  ## Heureka !
##
## 	Two-sample Kolmogorov-Smirnov test

## data:  sqrt(24000) * dip24k and sqrt(20000) * dip20k
## D = 0.0014, p-value = 0.292
## alternative hypothesis: two-sided
##

##--- Compute P-values of all pairwise KS tests: -- takes several minutes!!
m <- length(N.k.set)
P.vals <- matrix(NA, m,m, dimnames = list(names(N.k.set),names(N.k.set)))
diag(P.vals) <- 0
for(i in 1:(m-1)) {
    cat("i = ", i)
    d.i <- d.dip(N.k.set[i])
    for(j in (i+1):m) {
        cat(".")
        p <- ks.test(d.i, d.dip(N.k.set[j])) $ p.value
        P.vals[i,j] <- P.vals[j,i] <- p
    }; cat("\n")
}

round(P.vals, 5)
cP.vals <- format(round(P.vals,5))
cP.vals[upper.tri(P.vals)] <- ""
diag(cP.vals) <- "."
noquote(cP.vals)
##        8'000   12'000  16'000  20'000  24'000  32'000  36'000
## 8'000  .
## 12'000 0.00005 .
## 16'000 0.00000 0.01397 .
## 20'000 0.00000 0.00047 0.19526 .
## 24'000 0.00000 0.00006 0.03016 0.29203 .
## 32'000 0.00000 0.00000 0.00021 0.00034 0.02567 .
## 36'000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00154 .
##
## Hmm, not so good .. why is 36'000 different from 32'000 ??
save(P.vals, fd, fd., fdg, fdw,
     file = "asymp-res.rda")

(xl. <- extendrange(d.dip(20)))
xi <- seq(xl.[1], xl.[2], length.out = 2000)

plot(xi, Fn20k(xi) - Fn12k(xi), type = "l", col = 2)
## looks like the difference is largest where the density is large:
den20k <- density(sqrt(20000)* dip20k,
                  n = 2048, from = xl.[1], to = xl.[2], cut = 1)
f.scale <- .99 * par("usr")[3] / max(den20k$y)
with(den20k, lines(f.scale * y ~ x, col = "pink"))

## Graphical "KS test"
plot(xi, Fn24k(xi) - Fn20k(xi), type = "l", col = 2)
## looks much less systematic ...
f.scale <- .99 * par("usr")[3] / max(den20k$y)
with(den20k, lines(f.scale * y ~ x, col = "pink"))

## 36k  vs  32k  which looks "bad" in  KS test:
plot(xi, Fn36k(xi) - Fn32k(xi), type = "l", col = 2)
## looks much less systematic ...
f.scale <- .99 * par("usr")[3] / max(den20k$y)
with(den20k, lines(f.scale * y ~ x, col = "pink"))

##--- Ok, look even closer for systematic:
require(RColorBrewer)

Nmax <- max(N.k.set)# too large now
## rather just
Nmax <- 44
Fn <- ecdf(d.dip(Nmax))
Fn.xi <- Fn(xi)
ii <- seq_along(Ns <- N.k.set[N.k.set < Nmax])
opal <- palette(brewer.pal(length(Ns), "Set3"))# "Dark2" if  nColor <= 8

yrng <- range(Fn.xi - ecdf(d.dip(min(Ns)))(xi),
              Fn.xi - ecdf(d.dip(max(Ns)))(xi))
plot(range(xi), yrng, type = "n",
     ylab = "", xlab = "D_N  =^=  dip( runif(N) )",
     main = paste("F_n(D_{N=",Nmax,"'000})  -  F_n(D_N) ;   n=",
                  mFormat(length(d.dip(Nmax))), sep=''))
abline(h=0, col="gray")
for(i in ii)
    lines(xi, Fn(xi) - ecdf(d.dip(Ns[i]))(xi), col = i+1)
legend("right", legend=rev(paste("N = ", format(Ns), "'000", sep='')),
       col= rev(ii+1), lty=1, lwd=2, inset = .05)

palette(opal)

## How does the upper tail look?
## log(1 - P) = log(P{X >= x})
