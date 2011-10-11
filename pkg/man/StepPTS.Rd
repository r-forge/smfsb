\name{StepPTS}
\alias{StepPTS}
\title{Create a function for advancing the state of an SPN by using a simple approximate Poisson time stepping method}
\description{
This function creates a function for advancing the state of an SPN model using a simple approximate Poisson time stepping method. The resulting function (closure) can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of SPN models.
}
\usage{
StepPTS(N,dt=0.01)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic
    Petri net. Should contain \code{N$Pre}, a matrix representing the
    LHS stoichiometries, \code{N$Post}, a matrix representing the RHS
    stoichiometries, and \code{N$h}, a function representing the rates
    of the reaction processes. \code{N$h} should have first argument
    \code{x}, a vector representing the current state of the system, second argument
    \code{t}, a scalar representing the current simulation time (in the
    typical time-homogeneous case, \code{N$h} will ignore this argument).
    \code{N$h} may posses additional arguments, representing reaction rates, for example. \code{N} does not need to contain an initial marking, \code{N$M}. \code{N$M} will be ignored by most functions which use the resulting function closure.}
  \item{dt}{Time step to be used by the Poisson time stepping integration method. Defaults to 0.01.}
}
\value{
An R function which can be used to advance the state of the SPN model \code{N} by using a Poisson time stepping method with step size \code{dt}. The function closure has interface \code{function(x0,t0,deltat,...)}, where \code{x0} and \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepGillespie}}, \code{\link{StepCLE}}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# load up the LV model
data(spnModels)
# create a stepping function
stepLV=StepPTS(LV)
# step the function
print(stepLV(c(x1=50,x2=100),0,1))
# integrate the process and plot it
out = simTs(c(x1=50,x2=100),0,20,0.1,stepLV)
plot(out)
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
