\name{simSample}
\alias{simSample}
\title{Simulate a many realisations of a model at a given fixed time in the future given an initial time and state, using a function (closure) for advancing the state of the model}
\description{
This function 
simulates 
many realisations of a model at a given fixed time in the future given an initial time and state, using a function (closure) for advancing the state of the model
, such as created by \code{\link{StepGillespie}} or \code{\link{StepSDE}}.
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
\seealso{\code{\link{StepSDE}}, \code{\link{StepGillespie}}, \code{\link{simTimes}}, \code{\link{simTs}}}
\examples{
out3 = simSample(100,c(x1=50,x2=100),0,20,stepLVc)
hist(out3[,"x2"])
}

\keyword{smfsb}
