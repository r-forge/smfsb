# R wrapper for C code

stepLVc<-function(tt,x,deltat,th=c(1,0.005,0.6))
{
	new=.C("stepLV",as.double(tt),as.integer(x),as.double(deltat),as.double(th))[[2]]
	names(new)=names(x)
	new
}

# eof

