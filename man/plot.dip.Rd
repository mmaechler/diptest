\name{plot.dip}
\alias{plot.dip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.dip(x, do.points = (n < 20), colG = "red3", colL = "blue3", colM = "forest green", col.points = par("col"), col.hor = col.points, doModal = TRUE, doLegend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{do.points}{
%%     ~~Describe \code{do.points} here~~
}
  \item{colG}{
%%     ~~Describe \code{colG} here~~
}
  \item{colL}{
%%     ~~Describe \code{colL} here~~
}
  \item{colM}{
%%     ~~Describe \code{colM} here~~
}
  \item{col.points}{
%%     ~~Describe \code{col.points} here~~
}
  \item{col.hor}{
%%     ~~Describe \code{col.hor} here~~
}
  \item{doModal}{
%%     ~~Describe \code{doModal} here~~
}
  \item{doLegend}{
%%     ~~Describe \code{doLegend} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, do.points = (n < 20), colG = "red3", colL = "blue3", 
    colM = "forest green", col.points = par("col"), col.hor = col.points, 
    doModal = TRUE, doLegend = TRUE, ...) 
{
    stopifnot(is.integer(n <- x$n), is.numeric(D <- x$dip), length(lh <- x$lo.hi) == 
        2)
    Fn <- ecdf(x$x)
    cl <- x$call[1:2]
    cl[[1]] <- as.name("ecdf")
    names(cl)[2] <- ""
    attr(Fn, "call") <- cl
    chD <- formatC(D, digits = pmax(3, getOption("digits") - 
        2))
    tit <- bquote("Dip" ~ ~{
        D[n] == D[.(n)]
    } == .(chD))
    plot(Fn, do.points = do.points, col.points = col.points, 
        col.hor = col.hor, verticals = TRUE, col.vert = "sky blue", 
        lwd = 2, ...)
    title(tit, adj = 0, line = 1.25)
    aLine(x, "gcm", col = colG)
    aLine(x, "lcm", col = colL)
    if (doModal) {
        x12 <- x$x[lh]
        abline(v = x12, col = colM, lty = 2)
        op <- par(mgp = c(3, 1/16, 0))
        axis(3, at = x12, label = expression(x[L], x[U]), tick = FALSE, 
            col.axis = colM)
        par(op)
    }
    if (doLegend) 
        legend("topleft", bty = "n", c("greatest convex minorant GCM", 
            "least concave majorant    LCM"), col = c(colG, colL), 
            lty = 1, lwd = 1.5)
    invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
