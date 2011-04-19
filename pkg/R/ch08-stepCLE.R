# function

StepCLE=function(N,dt=0.01)
{
        S=t(N$Post-N$Pre)
	v=ncol(S)
	sdt=sqrt(dt)
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	       			h=N$h(tt,x,...)
				dw=rnorm(v,0,sdt)
				dx=S %*% (h*dt + sqrt(h)*dw)
       	         		x=x+as.vector(dx)
				tt=tt+dt
				if (tt > termt)
					return(x)
	       	 	}
		}
	)		
}	



# eof

