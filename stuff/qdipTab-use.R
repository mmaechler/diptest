####--- Use the simulated quantiles to compute P-values

### Now only working with qdipTab
library(diptest)
data(qDiptab)
(dnq <- dimnames(qDiptab))
(P.p <- as.numeric(dnq$Pr))## 0  0.01 ..... 0.99999 1
(n <- as.integer(dnq$n))   ## 4 5..10 15 ... 1000 2000 5000
## sqrt(n) * D  has limit distribution (sometimes..)
rq <- sqrt(n) * qDiptab

filled.contour(n, P.p, log10(rq),
               plot.title= contour(n, P.p, log10(rq), add = TRUE))
contour(n, P.p, logit(rq))

## correct P = {0,1} to {1/(2N), 1 - 1/(2N)} for N = 1000001 - 1 = 1e6
N <- 1000000
Pcp <- P.p
Pcp[P.p == 0] <-     1/(2*N)
Pcp[P.p == 1] <- 1 - 1/(2*N)

str(dd <- cbind(expand.grid(n=n, p=Pcp), rq = c(rq)))

coplot(rq ~ p | log10(n), data = dd)
coplot(p ~ rq | log10(n), data = dd)

## NOTA BENE:  qlogis(.) === logit(.)
coplot(qlogis(p) ~ rq | log10(n), data = dd)

library(mgcv)
g1 <- gam(qlogis(p) ~ s(log(n),rq), data = dd)
summary(g1) # 21.56 d.f. // Deviance explained 88.5%

if(FALSE) ## Not quite --- 50 warnings --- need to use  weights! etc
g2 <- gam(p ~ s(log(n),rq), family = binomial, data = dd)
