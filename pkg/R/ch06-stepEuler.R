# function

StepEuler=function(RHSfun,dt=0.01)
{
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	         		x=x+RHSfun(tt,x,...)*dt
				tt=tt+dt
				if (tt > termt)
					return(x)
	       	 	}
		}
	)		
}	



# eof

