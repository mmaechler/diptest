\name{qDiptab}
\alias{qDiptab}
\title{Table of Quantiles from a Large Simulation for Hartigan's Dip Test}
\docType{data}
\description{
  Whereas Hartigan(1985) published a table of empirical percentage
  points of the dip statistic (see \code{\link{dip}}) based on N=9999
  samples of size \eqn{n} from \eqn{U[0,1]}, our table of empirical
  quantiles is currently based on N=100'001 samples for each \eqn{n}.
}
\note{
  Taking N=100'001 ensures that all the \code{\link{quantile}(X, p)}
  used here are exactly order statistics \code{sort(X)[k]}.
}
\format{
  A numeric matrix %may change: of dimension 17 * 26,
  where each row corresponds to sample size \eqn{n}, and each column to
  a probability (percentage) in \eqn{[0,1]}.  The dimnames are coercable
  to these values, see the examples.
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}
}
\examples{
data(qDiptab)
str(qDiptab)
## the sample sizes `n' :
(nn <- as.integer(rownames(qDiptab)))
## the probabilities:
P.p <- as.numeric(print(colnames(qDiptab)))
nP <- length(P.p)
PdipRn <- qDiptab*sqrt(nn)
colP <- c(1:5,6,5:2,8:9,1)
matplot(nn, PdipRn, type="b", log="x", col = colP,
        main = "N^.5 * P-values for the 'dip' (H0=uniform)    vs. N")
mtext("(log scale for N)")
legend(2, 1.4,
       legend=rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
	 dimnames(qDiptab)[[2]], sep=": ")),
       col = rev(colP), lty = rev(rep(1:5,length=nP)),
       ncol = 4, bty="n", cex = 0.8)
}
\keyword{datasets}

