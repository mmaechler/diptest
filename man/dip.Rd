\name{dip}
\alias{dip}
\title{Compute Hartigan's Dip Test Statistic for Unimodality}
\description{
  Computes Hartigan's dip test statistic for testing unimodality,
  and additionally the modal interval.
}
\usage{
dip(x, full.result=FALSE, debug=FALSE)
}
\arguments{
  \item{x}{numeric; the data.}
  \item{full.result}{logical; if \code{TRUE} returns the full result
    list, see below.}
  \item{debug}{logical; if true, some tracing information is printed
    (from the C routine).}
}
\value{
  depending on \code{full.result} either, a number, the dip statistic or
  a list with components
  \item{dip}{the dip statistic}
  \item{xl}{lower end of modal interval}
  \item{xu}{upper end of modal interval}

  ...

  ...
}
\references{
  P. M. Hartigan (1985)
  Computation of the Dip Statistic to Test for Unimodality;
  \emph{Applied Statistics} \bold{34}, 320--325.

  J. A. Hartigan and P. M. Hartigan (1985)
  The Dip Test of Unimodality;
  \emph{Annals of Statistics} \bold{13}, 70--84.
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}, based on
  earlier code from Dario Ringach \email{dario@wotan.cns.nyu.edu}
}
\seealso{\code{\link{isoreg}} for isotonic regression.}
\examples{
data(exHartigan)
plot(density(exHartigan))
dip(exHartigan)

str(dip(exHartigan, full = TRUE, debug = TRUE))

data(faithful)
fE <- faithful$eruptions
plot(density(fE))
str(dip(fE, full = TRUE, debug = TRUE))

data(precip)
plot(density(precip))
str(dip(precip, full = TRUE, debug = TRUE))
}
\keyword{htest}
\keyword{distribution}
