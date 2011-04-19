\name{simTs}
\alias{simTs}
\title{Simulate a model on a regular grid of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates a model on a regular grid of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.
}
\usage{
simTs(t0=0,x0,tt=100,dt=0.1,stepFun,...)
}
\arguments{
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{x0}{The initial state of the process at time \code{t0}.}
  \item{tt}{The terminal time of the simulation.}
  \item{dt}{The time step of the output. Note that this time step relates only to the volume of the output, and has no bearing on the accuracy of the simulation process.}
  \item{stepFun}{A function (closure) for advancing the state of the process, such as produced by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R \code{\link{ts}} object representing the simulated process.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepGillespie}}, \code{\link{simTimes}}, \code{\link{simSample}}}
\examples{
# example
# Build a SPN for the Lotka-Volterra system
N=list()
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(t,x,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# create a stepping function
stepLV=StepGillespie(N)

# simulate a realisation of the process and plot it
out=simTs(0,c(x1=50,x2=100),100,0.1,stepLV)
plot(out)
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
