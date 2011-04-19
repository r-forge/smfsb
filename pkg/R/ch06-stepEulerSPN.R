# function

StepEulerSPN=function(N,dt=0.01)
{
        S=t(N$Post-N$Pre)
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	       			h=N$h(tt,x,...)
       	         		x=x+as.vector(S %*% h)*dt
				tt=tt+dt
				if (tt > termt)
					return(x)
	       	 	}
		}
	)		
}	



# eof

