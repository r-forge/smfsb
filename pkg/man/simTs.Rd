\name{simTs}
\alias{simTs}
\title{Simulate a model on a regular grid of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates single realisation of a model on a regular grid of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.
}
\usage{
simTs(x0,t0=0,tt=100,dt=0.1,stepFun,...)
}
\arguments{
  \item{x0}{The initial state of the process at time \code{t0}.}
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{tt}{The terminal time of the simulation.}
  \item{dt}{The time step of the output. Note that this time step
    relates only to the recorded output, and has no bearing on the accuracy of the simulation process.}
  \item{stepFun}{A function (closure) for advancing the state of the process, such as produced by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R \code{\link{ts}} object representing the simulated process.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepGillespie}},
  \code{\link{StepSDE}}, \code{\link{simTimes}}, \code{\link{simSample}}, \code{\link{as.timedData}}}
\examples{
# load the LV model
data(spnModels)
# create a stepping function
stepLV = StepGillespie(LV)
# simulate a realisation of the process and plot it
out = simTs(c(x1=50,x2=100),0,100,0.1,stepLV)
plot(out)
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
