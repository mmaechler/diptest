\name{qDiptab}
\alias{qDiptab}
\docType{data}
\title{Table of Quantiles from a Large Simulation for Hartigan's Dip Test}
\description{
  Whereas Hartigan(1985) published a table of empirical percentage
  points of his dip statistic (see \code{\link{dip}}) based on N=9999
  samples of size \eqn{n} from \eqn{U[0,1]}, our tables is currently
  based on N=40'000 samples for each \eqn{n}.
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}
}
\examples{
data(qDiptab)
str(qDiptab)
nP <- ncol(qDiptab)
(nn <- as.integer(row.names(qDiptab)))
PdipRn <- qDiptab*sqrt(nn)

colP <- c(1:5,6,5:2,8:9,1)
matplot(nn, PdipRn, type="b", log="x", col = colP,
        main = "N^.5 * P-values for the 'dip' (H0=uniform)    vs. N")
mtext("(log scale for N)")
legend(.88*min(nn), par("usr")[4],
       legend=rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
	 dimnames(qDiptab)[[2]], sep=": ")),
       col = rev(colP), lty = rev(rep(1:5,length=nP)),
       ncol = 4, bty="n", cex = 0.8)
}
\keyword{data}
