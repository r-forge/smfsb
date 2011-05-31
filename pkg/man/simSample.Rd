\name{simSample}
\alias{simSample}
\title{Simulate a many realisations of a model at a given fixed time in the future given an initial time and state, using a function (closure) for advancing the state of the model}
\description{
This function 
simulates 
many realisations of a model at a given fixed time in the future given an initial time and state, using a function (closure) for advancing the state of the model
, such as created by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.
}
\usage{
simSample(n=100,x0,t0=0,deltat,stepFun,...)
}
\arguments{
  \item{n}{The number of samples required. Defaults to 100.}
  \item{x0}{The initial state of the process at time \code{t0}.}
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{deltat}{The amount of time in the future of \code{t0} at which samples of the system state are required.}
  \item{stepFun}{A function (closure) for advancing the state of the process, such as produced by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R matrix whose rows represent the simulated states of the process at time \code{t0+deltat}.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepGillespie}}, \code{\link{simTimes}}, \code{\link{simTs}}}
\examples{
# example
# Build a SPN for the Lotka-Volterra system
N=list()
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(x,t,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# create a stepping function
stepLV = StepGillespie(N)

# simulate a sample from the transition kernel
out3 = simSample(100,c(x1=50,x2=100),0,20,stepLV)
hist(out3[,"x2"])
}

\keyword{smfsb}
