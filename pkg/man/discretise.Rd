\name{discretise}
\alias{discretise}
\title{Discretise output from a discrete event simulation algorithm}
\description{
This function discretise output from a discrete event simulation algorithm such as the Gillespie algorithm onto a regular time grid, and returns the results as an R ts object.
}
\usage{
discretise(out, dt=1, start=0)
}
\arguments{
  \item{out}{A list containing discrete event simulation output in the form of that produced by \code{\link{gillespie}}.}
  \item{dt}{The time step required for the output of the discretisation process. Defaults to one time unit.}
  \item{start}{The start time for the output. Defaults to zero.}
}
\value{
An R \code{\link{ts}} object containing the discretised output.
}
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}}, \code{\link{gillespie}}, \code{\link{gillespied}}, \code{\link{ts}}}
\examples{
# example

# Build a SPN for the Lotka-Volterra system
N=list()
N$M=c(x1=50,x2=100)
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(t,x,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# simulate a realisation of the process and plot it
out=gillespie(N,10000)
op=par(mfrow=c(2,2))
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")

# use the "discretise" function to map it to an R "ts" object
plot(discretise(out,dt=0.01),plot.type="single",lty=c(1,2))
par(op)
}

\keyword{smfsb}
