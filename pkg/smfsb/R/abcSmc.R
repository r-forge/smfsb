## abcSmc.R


require(parallel)

abcSmc <- function(N, rprior, dprior, rdist, rperturb, dperturb, factor=10,
                   steps=15, verb=FALSE) {
    priorLW = log(rep(1/N, N))
    priorSample = mclapply(as.list(priorLW), function(x) {rprior()})
    for (i in steps:1) {
        if (verb) message(paste(i,""), appendLF=FALSE)
        out = abcSmcStep(dprior, priorSample, priorLW, rdist, rperturb,
                         dperturb, factor)
        priorSample = out[[1]]
        priorLW = out[[2]]
    }
    if (verb) message("")
    t(sapply(sample(priorSample, N, replace=TRUE, prob=exp(priorLW)), identity))
}

abcSmcStep <- function(dprior, priorSample, priorLW, rdist, rperturb,
                       dperturb, factor=10) {
    n = length(priorSample)
    mx = max(priorLW)
    rw = exp(priorLW - mx)
    prior = sample(priorSample, n*factor, replace=TRUE, prob=rw)
    prop = mcMap(rperturb, prior)
    dist = mcMap(rdist, prop) # forward simulation in parallel
    qCut = quantile(unlist(dist), 1/factor)
    new = prop[dist < qCut]
    lw = mcMap( function(th) {
        terms = priorLW + sapply(priorSample,
                                 function(x){dperturb(th, x, log=TRUE)})
        mt = max(terms)
        denom = mt + log(sum(exp(terms - mt)))
        dprior(th, log=TRUE) - denom
    } , new)
    lw = unlist(lw)
    mx = max(lw)
    lw = lw - mx
    nlw = log(exp(lw)/sum(exp(lw)))
    list(sample = new, lw = nlw)
}


## eof

