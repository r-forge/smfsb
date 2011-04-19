# function

StepGillespie=function(N)
{
        S=t(N$Post-N$Pre)
        v=ncol(S)
	return(
		function(tt,x,deltat,...)
		{
			termt=tt+deltat
	        	repeat {
       	       			h=N$h(tt,x,...)
				h0=sum(h)
		                if (h0<1e-10)
               			        tt=1e99
                		else
                        		tt=tt+rexp(1,h0)
				if (tt>=termt)
					return(x)
       	       			j=sample(v,1,prob=h)
       	         		x=x+S[,j]
	       	 	}
		}
	)		
}	

simTs=function(t0=0,x0,tt=100,dt=0.1,stepFun,...)
{
	n=(tt-t0) %/% dt + 1
	u=length(x0)
	names=names(x0)
	mat=matrix(nrow=n,ncol=u)
	x=x0
	t=t0
	mat[1,]=x
	for (i in 2:n) {
		t=t+dt
		x=stepFun(t,x,dt,...)
		mat[i,]=x
	}
	ts(mat,start=t0,deltat=dt,names=names)
}

simTimes=function(t0=0,x0,times,stepFun,...)
{
	n=length(times)
	u=length(x0)
	names=names(x0)
	mat=matrix(nrow=n,ncol=u)
	x=x0
	t=t0
	for (i in 1:n) {
		x=stepFun(t,x,times[i]-t,...)
		t=times[i]
		mat[i,]=x
	}
	rownames(mat)=times
	colnames(mat)=names
	mat
}

simSample=function(n=100,t0=0,x0,deltat,stepFun,...)
{
	u=length(x0)
	names=names(x0)
	mat=matrix(nrow=n,ncol=u)
	for (i in 1:N) {
		mat[i,]=stepFun(t0,x0,deltat,...)
	}
	colnames(mat)=names
	mat
}

# example

# eof

