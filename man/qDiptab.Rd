\name{P.dip}
\alias{P.dip}
\title{P Value Table from a Large Simulation for Hartigan's Dip Test}
\description{

}
\examples{
data(P.dip)
str(P.dip)
nP <- ncol(P.dip)
(nn <- as.integer(row.names(P.dip)))
P.dip.Rn <- P.dip*sqrt(nn) 

colP <- c(1:5,6,5:2,8:9,1)
matplot(nn, P.dip.Rn, type="b", log="x", col = colP,
        main = "N^.5 * P-values for the 'dip' (H0=uniform)    vs. N")
mtext("(log scale for N)")
legend(.88*min(nn), par("usr")[4], 
       legend=rev(paste(c(paste(c(1:9,0)),letters)[1:nP],
	 dimnames(P.dip)[[2]], sep=": ")),
       col = rev(colP), lty = rev(rep(1:5,length=nP)),
       ncol = 4, bty="n", cex = 0.8)
}
\keyword{data}
