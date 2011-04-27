# function

require(deSolve)

StepODE=function(RHSfun)
{
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
			res=deSolve::ode(x,times=c(tt,termt),func=function(t,x,...) list(RHSfun(t,x,...)),...)
			res=res[2,]
			res[2:length(res)]
		}
	)		
}	



# eof

