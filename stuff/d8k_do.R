stopifnot(require("diptest"))
setwd("/u/maechler/R/Pkgs/diptest/stuff")

N <- 8000
n.sim <- 1000001
d8k <- numeric(n.sim)
.pt <- proc.time()

set.seed(123)
for(i in 1:n.sim) {
    d8k[i] <- dip(runif(N))
    if(i %% 100 == 0) {
        cat(".")
        if(i %% 10000 == 0) {
            cat("",i,"\n")
            .opt <- .pt ; .pt <- proc.time(); .pt - .opt
        }
    }
}
