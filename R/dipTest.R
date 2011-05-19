dip.test <- function(x, simulate.p.value = FALSE, B = 2000) {
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    stopifnot(is.numeric(x))
    n <- length(x) # *is* integer
    D <- dip(x)
    if(simulate.p.value) {
        P <- mean(D <= replicate(B, dip(runif(n))))
    } else {
        data(qDiptab, package="diptest")
        dn <- dimnames(qDiptab)
        nn <- as.integer(dn[["n"]])
        P.s <- as.numeric(noquote(dn[["Pr"]]))
}
    structure(list(statistic = c(D = D),
                   p.value = sw$pw,
                   method = "Hartigans' test for unimodality",
                   data.name = DNAME)
              class = "htest")
}
