# function

StepFRM=function(N)
{
        S=t(N$Post-N$Pre)
        v=ncol(S)
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	       			h=N$h(tt,x,...)
				pu=rexp(v,h)
				j=which.min(pu)
				tt=tt+pu[j]
				if (tt>=termt)
					return(x)
       	         		x=x+S[,j]
	       	 	}
		}
	)		
}	


# eof

