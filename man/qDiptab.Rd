\name{qDiptab}
\alias{qDiptab}
\title{Table of Quantiles from a Large Simulation for Hartigan's Dip Test}
\docType{data}
\description{
  Whereas Hartigan(1985) published a table of empirical percentage
  points of the dip statistic (see \code{\link{dip}}) based on N=9999
  samples of size \eqn{n} from \eqn{U[0,1]}, our table of empirical
  quantiles is currently based on N=1'000'001 samples for each \eqn{n}.
}
\note{
  Taking N=1'000'001 ensures that all the \code{\link{quantile}(X, p)}
  used here are exactly order statistics \code{sort(X)[k]}.
}
\format{
  A numeric matrix %may change: of dimension 17 * 26,
  where each row corresponds to sample size \eqn{n}, and each column to
  a probability (percentage) in \eqn{[0,1]}.  The dimnames are coercable
  to these values, see the examples.
}
\seealso{\code{\link{dip}}, also for the references.}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}
}
\examples{
data(qDiptab)
str(qDiptab)
## the sample sizes `n' :
(nn <- as.integer(rownames(qDiptab)))
## the probabilities:
P.p <- as.numeric(print(colnames(qDiptab)))

## This is as "Table 1" in Hartigan & Hartigan (1985) -- but more accurate
ps <- c(1,5,10,50,90,95,99, 99.5, 99.9)/100
tab1 <- qDiptab[nn <= 200,  as.character(ps)]
round(tab1, 4)
}
\keyword{datasets}
