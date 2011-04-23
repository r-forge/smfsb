# R wrapper for C code

stepLVc<-function(t0,x0,dt,th=c(1,0.005,0.6))
{
	new=.C("stepLV",as.double(t0),as.integer(x0),as.double(dt),as.double(th),PACKAGE="smfsb")[[2]]
	names(new)=names(x0)
	new
}

# eof

