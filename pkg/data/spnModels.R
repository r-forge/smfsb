# Example SPN Models

# SPN for the Lotka-Volterra system
LV=list()
LV$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
LV$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
LV$M=c(x1=50,x2=100)
LV$h=function(x,t,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# Imigration-Death process
ID=list()
ID$Pre=matrix(c(0,1),ncol=1)
ID$Post=matrix(c(1,0),ncol=1)
ID$h=function(x,t,th=c(lambda=1,mu=0.1))
{
 with(as.list(c(x,th)),{
  return(c(lambda,mu*x))
 })
}

# Birth-Death process
BD=list()
BD$Pre=matrix(c(1,1),ncol=1)
BD$Post=matrix(c(2,0),ncol=1)
BD$h=function(x,t,th=c(lambda=1,mu=1.1))
{
 with(as.list(c(x,th)),{
  return(c(lambda*x,mu*x))
 })
}

Dimer=list()
Dimer$Pre=matrix(c(2,0,0,1),ncol=2,byrow=TRUE)
Dimer$Post=matrix(c(0,1,2,0),ncol=2,byrow=TRUE)
Dimer$h=function(x,t,th=c(th1=1.66e-3,th2=0.2)) {
 with(as.list(c(x,th)),{
        return(c(th1*x1*(x1-1)/2, th2*x2))
	})	
    }




# eof 

