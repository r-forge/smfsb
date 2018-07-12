## abcRun.R

require(parallel)

abcRun <- function(n, rprior, rdist) {
    v = vector("list", n)
    p = mclapply(v, function(x){ rprior() }) # use mcMap instead??
    d = mclapply(p, rdist) # forward simulation in parallel
    pm = t(sapply(p, identity))
    if (dim(pm)[1] == 1) pm = as.vector(pm)
    dm = t(sapply(d, identity))
    if (dim(dm)[1] == 1) dm = as.vector(dm)
    list(param=pm, dist=dm)
}

## eof

