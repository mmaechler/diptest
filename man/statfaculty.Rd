\name{statfaculty}
\alias{statfaculty}
\title{Faculty Quality in Statistics Departments}
\docType{data}
\usage{
data(statfaculty)
}
\description{
 Faculty quality in statistics departments was assessed as part of a
 larger study reported by Scully(1982).
}
\format{A numeric vector of 63 (integer) numbers, sorted increasingly,
  as reported by the reference.
}
\references{
  J. A. Hartigan and P. M. Hartigan (1985)
  The Dip Test of Unimodality;
  \emph{Annals of Statistics} \bold{13}, 70--84.
}
\source{
  M. G. Scully (1982)
  Evaluation of 596 programs in mathematics and physical sciences;
  \emph{Chronicle Higher Educ.} \bold{25} 5, 8--10.
}
\examples{
data(statfaculty)
plot(dH <- density(statfaculty))
rug(jitter(statfaculty))
}
\keyword{datasets}
