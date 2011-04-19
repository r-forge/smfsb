\name{simTimes}
\alias{simTimes}
\title{Simulate a model at a specified set of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates a model at a specified set of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.
}
\usage{
simTimes(t0=0,x0,times,stepFun,...)
}
\arguments{
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{x0}{The initial state of the process at time \code{t0}.}
  \item{times}{A vector of times at which the state of the process is required. It is assumed that the times are in increasing order, and that the first time is at least as big as \code{t0}.}
  \item{stepFun}{A function (closure) for advancing the state of the process, such as produced by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R matrix where each row represents the state of the process at one of the required times.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepGillespie}}, \code{\link{simTs}}, \code{\link{simSample}}}
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

# simulate a realisation using simTimes
times=seq(0,100,by=0.1)
plot(ts(simTimes(0,c(x1=50,x2=100),times,stepLV),start=0,deltat=0.1),plot.type="single",lty=1:2)
# simulate a realisation at irregular times
times=c(0,10,20,50,100)
out2=simTimes(0,c(x1=50,x2=100),times,stepLV)
print(out2)
}

\keyword{smfsb}
