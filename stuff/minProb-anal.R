#### Analysis of the result of ./sim-minProb.R
#### ~~~~~~~~                  ~~~~~~~~~~~~~~~

setwd("/u/maechler/R/Pkgs/diptest/stuff")
(nMin <- readRDS(file= "minProb.rds"))
Ns <- attr(nMin, "Ns")
nn <- as.integer(names(nMin))
nMin[nMin > 0] ## now only goes up to n=13,  even for Ns == 500'000

## Actually  nMin ~ Bin(pm(n), Ns) --> estimate pm(n) parametrically:
##           =====================             ======
glm.n <- glm(cbind(nMin,Ns) ~ log(nn), family = binomial)
summary(glm.n)
## ^^ shows that logistic link  &  log(nn)  was not good enough:

glm.s <- glm(cbind(nMin,Ns) ~ sqrt(nn), family = binomial)
summary(glm.s)
glm.is <- glm(cbind(nMin,Ns) ~ I(nn^-0.5), family = binomial) #  "best"
summary(glm.is)
glm.I <- glm(cbind(nMin,Ns) ~ I(nn^-1), family = binomial)
summary(glm.I)

matplot(nn, cbind(nMin
                  ,Ns * predict(glm.s, type="response") # "2"
                  ,Ns * predict(glm.n, type="response") # "3"
                  ,Ns * predict(glm.is, type="response")# "4" (best!)
                  ,Ns * predict(glm.I, type="response") # "5"
                  ), type='b',
        log = "xy",
        ylab = "nMin  & 4 different glm() predictions")

## actually with large simulation, "glm.is" are not bad
plot(nMin ~ nn, type = 'o')

## or
require(gplots)
plotCI(nn, nMin, uiw = sqrt(nMin*(Ns-nMin)/Ns), type = 'o')

## in Log-Log scale  -- + error bars
plotCI(nn, nMin, uiw = 1.96*sqrt(nMin*(Ns-nMin)/Ns), type = 'o', log = 'xy')
## 2011-05-14: AAARGH the above now gives an error!
## whereas this still works: "was"
plot(nMin ~ nn, type = 'o', log = 'xy')
axis(2, at=c(2000,5000,10000))
axis(1, at=nn)

## Better confidence intervals (this takes time!):
ci <- Ns*sapply(nMin, function(k)
                binom.test(k, Ns, conf.level = .99)$conf.int)

plot(nMin ~ nn, type = 'o', log = 'xy', xlab = 'n', cex = 0.8)
segments(nn, ci[1,], nn, ci[2,], col = "tomato")
axis(2, at=c(2000,5000,10000))
axis(1, at=nn)
title(paste("#{dip = 1/(2 n)}  for Ns=",formatC(Ns,form="d"),"  samples"))
mtext("exact binomial 99%-confidence intervals", col = "tomato")
text(nn, nMin, formatC(signif(nMin/Ns, 2)), adj = -0.2, cex = 0.8)

## i.e., for small 'n'  the probability is remarkably high,
## where as for n = 5000,
## it is only visible if you look very carefully :

##___________ No longer ____ after fixing the 2003-10  bug ______

### for  n=5000 "5k":
load("dip5k.rda")# d5k has 100'000 values dip(runif(5000))

##     main = substitute("density of  " * P(D[n] == frac(1,2*n))

plot(density(d5k),
     main = substitute("density of  " * {
         D[n] == "dip(runif(5000)),   simulated  " } *
     N == nnn *" times", list(nnn = length(d5k))))
## or with `sfsmisc' package:
tkdensity(d5k, do.rug=FALSE, from.f = -4)
## --> Title of a paper: The distribution of the dip is "bimodal"  !!
##
## or  "The Dip Distribution has a Dip" !

