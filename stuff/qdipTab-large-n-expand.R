Ns <- 1000001

stopifnot(require("diptest") && is.character(data(qDiptab)))

## begin { from  ./asymp-distrib.R } -------------------------------------------
## These all have  n.sim = 1000'001  samples of
## dip(runif(N)) for "large N"
## where produced by scripts such as  ./d20k_do.R
N.k.set <- c(8,12,16,20,24,32,36,40)
## or automatically
patt <- "^dip(.*)k\\.rda"
N.k.set <- sort(as.integer(sub(patt, "\\1", list.files(pattern = patt))))
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
##----end{ from ./asymp-distrib.R } -------------------------------------------

##
dn <- dimnames(qDiptab)
nn  <- as.numeric(dn[[1]])
P.p <- as.numeric(dn[[2]])
names(N.k.set) <- 1000* N.k.set
qNew <- sapply(N.k.set, function(.)
               quantile(d.dip(., scaleUp=FALSE),
                        probs = P.p))
signif(qNew, 3)
qDiptab.n <- rbind(qDiptab,
                   t(qNew[, as.character(1000*c(10,20,40,72))]))
## the dimnames-names got lost in rbind():
names(dimnames(qDiptab.n)) <- names(dimnames(qDiptab))

attr(qDiptab.n, "N_1") <- as.integer(Ns - 1)
## --- Here comes the *NEW* one:
attach("~/R/Pkgs/diptest/data/qDiptab.rda")

qDiptab.prev <- qDiptab
qDiptab <- qDiptab.n
if(FALSE) # do not do this "accidentally" !
save(qDiptab, file="~/R/Pkgs/diptest/data/qDiptab.rda")
