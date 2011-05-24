# function

StepSDE <- function (drift, diffusion, dt = 0.01) 
{
    sdt = sqrt(dt)
    return(function(tt, x, deltat, ...) {
        termt = tt + deltat
	v = length(x)
        repeat {
	    dw = rnorm(v,0,sdt)
            x = x + drift(tt, x, ...)*dt + as.vector(diffusion(tt,x,...)%*%dw)
            tt = tt + dt
            if (tt > termt)
		return(x)
        }
    })
}


# eof

