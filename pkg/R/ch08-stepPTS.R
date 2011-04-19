# function

StepPTS=function(N,dt=0.01)
{
        S=t(N$Post-N$Pre)
	v=ncol(S)
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	       			h=N$h(tt,x,...)
				r=rpois(v,h*dt)
       	         		x=x+as.vector(S %*% r)
				tt=tt+dt
				if (tt > termt)
					return(x)
	       	 	}
		}
	)		
}	



# eof

