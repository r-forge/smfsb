\name{simTimes}
\alias{simTimes}
\title{Simulate a model at a specified set of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates a single realisation from a Markovian model and
records the state at a specified set of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.
}
\usage{
simTimes(x0,t0=0,times,stepFun,...)
}
\arguments{
  \item{x0}{The initial state of the process at time \code{t0}.}
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{times}{A vector of times at which the state of the process is required. It is assumed that the times are in increasing order, and that the first time is at least as big as \code{t0}.}
  \item{stepFun}{A function (closure) for advancing the state of the process, such as produced by \code{\link{StepGillespie}} or \code{\link{StepEulerSPN}}.}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R matrix where each row represents the state of the process at one of the required times. The row names contain the sampled times.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepGillespie}},
  \code{\link{simTs}}, \code{\link{simSample}},
  \code{\link{as.timedData}}, \code{\link{pfMLLik}}}
\examples{
# load the LV model
data(spnModels)
# create a stepping function
stepLV = StepGillespie(LV)
# simulate a realisation using simTimes
times = seq(0,100,by=0.1)
plot(ts(simTimes(c(x1=50,x2=100),0,times,stepLV),start=0,deltat=0.1),plot.type="single",lty=1:2)
# simulate a realisation at irregular times
times = c(0,10,20,50,100)
out2 = simTimes(c(x1=50,x2=100),0,times,stepLV)
print(out2)
}

\keyword{smfsb}
