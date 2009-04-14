##-*- mode: R; kept-new-versions: 21; kept-old-versions: 12; -*-

#### Originally, from
####   http://www.stat.washington.edu/wxs/Stat593-s03/Code/jeremy-unimodality.R
####
####-----------------------------------------------------------------------------------
## Diagnostic plots for clustering and the DIP test for unimodality
## Code written by Jeremy Tantrum, Winter 2003
##=================================================================

#### --> see ./jeremy-unimodality.R

## An example of it working: Olive oil data - region 2 - are areas 5 and 6
## different.

## -- MM:  Aargh, 'olives' *differs* from R package to R package
data(olives, package="classifly")
str(olives.CF <- olives)
data(olives, package="TWIX")
str(olives.TW <- olives)
if(FALSE)## These are completely different:
    data(oliveoil, package="pls")


## First, the meaning of  'Area' and 'Region' is swapped between the two:
with(olives.TW, table(Region, Area))
##                  Area
## Region            North Sardinia South
##   Calabria            0        0    56
##   Coast Sardinia      0       33     0
##   East Liguria       50        0     0
##   Inland Sardinia     0       65     0
##   North Apulia        0        0    25
##   Sicily              0        0    36
##   South Apulia        0        0   206
##   Umbria             51        0     0
##   West Liguria       50        0     0

with(olives.CF, table(Area, Region))
##                  Region
## Area                1   2   3
##   Calabria         56   0   0
##   Coast-Sardinia    0  33   0
##   East-Liguria      0   0  50
##   Inland-Sardinia   0  65   0
##   North-Apulia     25   0   0
##   Sicily           36   0   0
##   South-Apulia    206   0   0
##   Umbria            0   0  51
##   West-Liguria      0   0  50


## Not at all the same... aargh
{ op <- par(ask=TRUE)
  for(n in names(oo)) { plot(oo[,n], olives[,n], main=n) }
  par(op) }

## But looking at subsets
subset(olives.TW, Region == "Coast Sardinia")[, 1:8]
subset(olives.CF, Area   == "Coast-Sardinia")[, -(1:2)]

### See ./Stuetzle-stat593-S2003-olive.doc :
olives.WS <-
   read.table("/u/maechler/R/Pkgs/diptest/stuff/Stuetzle-stat593-S2003-olive.tab",
              header=TRUE)

## together with code in ./jeremy-unimodality.R :
all.equal(unname(olives.CF[,-2]), unname(olives.WS[,-2])) ## TRUE
## but clearly, also the last three *Variable* names are confused...  [aargh]

## How can we *match* the two sets?  Look at the [,2] variables
with(olives.CF, table(Area,  Region))
with(olives.WS, table(Class, Region))

## -> manually :
cl2area <- c(5,1,7,6, 4,2, 3,9,8)

## but can do this programmatically:
(areas <- levels(olives.CF[,2]))
c2a <- unique(cbind(olives.CF[,2], # $ Area
                    olives.WS[,2]))# $ Class
(c2a <- c2a[order(c2a[,2]) , ])
c2a <- c2a[,1]
stopifnot(all(c2a == cl2area),
          all(olives.CF[,2] == areas[c2a[olives.WS[,2]]]))

## Now to work with Jeremy's code, assuming he used Stuetzle's version of the data:
olive.area   <- olives.WS $ Class
olive.region <- olives.WS $ Region
i.rest <- !(names(olives.WS) %in% c("Class", "Region"))
olive        <- data.matrix(olives.WS[, i.rest])

##
x.labs <- olive.area[olive.region==2]
library(MASS)
g <- lda(olive[olive.region==2,], x.labs)
g.proj <- unclass(g)$scaling
##
x <- olive[olive.region==2,] %*% g.proj
plot.ucdf(x)
plot.silverman(x)
str(x.dip <- calcdip(x, plot.it=FALSE)) #-> . $dip = 0.149
dips <- rep(0,100)
for(i in 1:100) {
  x.boot <- unisample(x.dip$unicurve,length(x))
  dips[i] <- calcdip(x.boot, plot.it=FALSE)$dip
}
(p.value <- sum(dips>x.dip$dip)/100) # 0  ( < 0.01 )
# < 0.01
##
