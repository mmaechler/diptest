\name{statfaculty}
\alias{statfaculty}
\title{Faculty Quality in Statistics Departments}
\description{
 Faculty quality in statistics departments was assessed as part of a
 larger study reported by Scully(1982).
 We use the 63 (integer) numbers reported by the first reference.
}
\references{
  J. A. Hartigan and P. M. Hartigan (1985)
  The Dip Test of Unimodality;
  \emph{Annals of Statistics} \bold{13}, 70--84.

  M. G. Scully (1982)
  Evaluation of 596 programs in mathematics and physical sciences;
  \emph{Chronicle Higher Educ.} \bold{25} 5, 8--10.
}
\examples{
data(statfaculty)
plot(dH <- density(statfaculty))
rug(jitter(statfaculty))
}
\keyword{data}
