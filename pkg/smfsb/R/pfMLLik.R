## pfMLLik
## particle filter unbiased estimate of marginal likelihood (logged)

pfMLLik = function (n, simx0, t0, stepFun, dataLik, data) 
{
    times = c(t0, as.numeric(rownames(data)))
    deltas = diff(times)
    return(function(...) {
        xmat = simx0(n, t0, ...)
        ll = 0
        for (i in 1:length(deltas)) {
            xmat[] = t(apply(xmat, 1, stepFun, t0 = times[i], deltat = deltas[i], 
                ...))
            lw = apply(xmat, 1, dataLik, t = times[i + 1], y = data[i, 
                ], log = TRUE, ...)
            m = max(lw)
            rw = lw - m
            sw = exp(rw)
            ll = ll + m + log(mean(sw))
            rows = sample(1:n, n, replace = TRUE, prob = sw)
            xmat[] = xmat[rows,]
        }
        ll
    })
}


## eof

