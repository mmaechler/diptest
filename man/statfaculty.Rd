\name{exHartigan}
\alias{exHartigan}
\title{Hartigan's Artificial n-modal Example Data Set}
\description{
 63 (integer) numbers; unimodal or bimodal, that's the question.
}
\examples{
data(exHartigan)
plot(dH <- density(exHartigan))
rug(jitter(exHartigan))
}
\keyword{data}
