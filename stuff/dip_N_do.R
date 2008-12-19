stopifnot(exists("N"),     is.numeric(N),     length(N) == 1, N == round(N),
          exists("n.sim"), is.numeric(n.sim), length(n.sim) == 1,
          n.sim == round(n.sim))
stopifnot(require("diptest"))
setwd("/u/maechler/R/Pkgs/diptest/stuff")

dd <- numeric(n.sim)
.pt <- proc.time()

set.seed(12)
for(i in 1:n.sim) {
    dd[i] <- dip(runif(N))
    if(i %% 100 == 0) {
        cat(".")
        if(i %% 10000 == 0) {
            cat("",i,"\n")
            .opt <- .pt ; .pt <- proc.time(); .pt - .opt
        }
    }
}

nam <- paste("dip", floor(N/1000),"k", sep='')
assign(nam, dd)
(outf <- paste(nam, "rda", sep="."))
save(list = c(nam, "N"), file = outf)
