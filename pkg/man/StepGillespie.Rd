\name{StepGillespie}
\alias{StepGillespie}
\title{Create a function for advancing the state of an SPN by using the Gillespie algorithm}
\description{
This function creates a function for advancing the state of an SPN model using the Gillespie algorithm. The resulting function (closure) can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of SPN models.
}
\usage{
StepGillespie(N)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic
    Petri net (SPN). Should contain \code{N$Pre}, a matrix representing
    the LHS stoichiometries, \code{N$Post}, a matrix representing the
    RHS stoichiometries, and \code{N$h}, a function representing the
    rates of the reaction processes. \code{N$h} should have
    first argument \code{x}, a vector representing the
    current state of the system, and
    second
    argument \code{t}, a scalar representing the current simulation time
    (in the typical time-homogeneous case, \code{N$h} will ignore this
    argument).
    \code{N$h} may possess additional arguments, representing reaction rates, for example. \code{N} does not need to contain an initial marking, \code{N$M}. \code{N$M} will be ignored by most functions which use the resulting function closure.}
}
\value{
An R function which can be used to advance the state of the SPN model \code{N} by using the Gillespie algorithm. The function closure has interface \code{function(x0,t0,deltat,...)}, where \code{x0} and \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{gillespied}},
  \code{\link{simTs}}, \code{\link{simTimes}}, \code{\link{simSample}}, \code{\link{StepFRM}},
  \code{\link{StepPTS}}, \code{\link{StepCLE}}}
\examples{
# load up the Lotka-Volterra (LV) model
data(spnModels)
LV
# create a stepping function
stepLV = StepGillespie(LV)
# step the function
print(stepLV(c(x1=50,x2=100),0,1))
# simulate a realisation of the process and plot it
out = simTs(c(x1=50,x2=100),0,100,0.1,stepLV)
plot(out)
plot(out,plot.type="single",lty=1:2)
# simulate a realisation using simTimes
times = seq(0,100,by=0.1)
plot(ts(simTimes(c(x1=50,x2=100),0,times,stepLV),start=0,deltat=0.1),plot.type="single",lty=1:2)
# simulate a realisation at irregular times
times = c(0,10,20,50,100)
out2 = simTimes(c(x1=50,x2=100),0,times,stepLV)
print(out2)
}

\keyword{smfsb}
