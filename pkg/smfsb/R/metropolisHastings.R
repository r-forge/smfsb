## metropolisHastings.R


metropolisHastings <- function(init, logLik, rprop, dprop=function(new, old, ...){1},
                               dprior=function(x, ...){1}, iters=10000, thin=10,
                               verb=TRUE, debug=FALSE) {
    p=length(init)
    ll=-Inf
    mat=matrix(0,nrow=iters,ncol=p)
    colnames(mat)=names(init)
    x=init
    if (verb) message(paste(iters,"iterations"))
    for (i in 1:iters) {
	if (verb) message(paste(i,""),appendLF=FALSE)
	for (j in 1:thin) {
            prop=rprop(x)
            if (dprior(prop, log=TRUE) > -Inf) {
                llprop=logLik(prop)
                a = (llprop - ll +
                     dprior(prop, log=TRUE) - dprior(x, log=TRUE) +
                     dprop(x, prop, log=TRUE) - dprop(prop, x, log=TRUE))
                if (debug) {
                    message(paste0("x=",x,", prop=",prop,collapse="\n"))
                    message(paste0("ll=",ll,", llprop=",llprop,", a=",a))
                }
                if (log(runif(1)) < a) {
                    x=prop
                    ll=llprop
                }
            }
        }
	mat[i,]=x
    }
    if (verb) message("Done.")
    mat
}


## eof


