\name{pval.dip}
\alias{pval.dip}
\docType{data}
\title{P Value Table from a Large Simulation for Hartigan's Dip Test}
\description{
  todo
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}
}
\examples{
data(pval.dip)
str(pval.dip)
nP <- ncol(pval.dip)
(nn <- as.integer(row.names(pval.dip)))
PdipRn <- pval.dip*sqrt(nn)

colP <- c(1:5,6,5:2,8:9,1)
matplot(nn, PdipRn, type="b", log="x", col = colP,
        main = "N^.5 * P-values for the 'dip' (H0=uniform)    vs. N")
mtext("(log scale for N)")
legend(.88*min(nn), par("usr")[4],
       legend=rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
	 dimnames(pval.dip)[[2]], sep=": ")),
       col = rev(colP), lty = rev(rep(1:5,length=nP)),
       ncol = 4, bty="n", cex = 0.8)
}
\keyword{data}
