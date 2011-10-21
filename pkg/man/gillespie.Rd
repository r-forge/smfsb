\name{gillespie}
\alias{gillespie}
\title{Simulate a sample path from a stochastic kinetic model described by a stochastic Petri net}
\description{
This function simulates a single realisation from a discrete stochastic
kinetic model described by a stochastic Petri net (SPN).
}
\usage{
gillespie(N, n, ...)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic
    Petri net (SPN). Should contain \code{N$M}, a vector representing
    the initial marking of the net, \code{N$Pre}, a matrix representing
    the LHS stoichiometries, \code{N$Post}, a matrix representing the
    RHS stoichiometries, and \code{N$h}, a function representing the
    rates of the reaction processes. \code{N$h} should have
    first argument \code{x}, a vector representing the current state of
    the system, and
    second argument \code{t}, a scalar representing the current
    simulation time (in the typical time-homogeneous case, \code{N$h}
    will ignore this argument).
    \code{N$h} may posses additional arguments, representing reaction rates, for example.}
  \item{n}{An integer representing the number of events to simulate, excluding the initial state, \code{N$M}.}
  \item{...}{Additional arguments (such as reactions rates) will be passed into the function \code{N$h}.}
}
\value{
A list with first component \code{t}, a vector of length \code{n} containing event times and second component \code{x}, a matrix with \code{n+1} rows containing the state of the system. The \code{i}th row of \code{x} contains the state of the system prior to the \code{i}th event.
}
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}},
  \code{\link{discretise}}, \code{\link{gillespied}}, \code{\link{StepGillespie}}}
\examples{
# load the LV model
data(spnModels)
# simulate a realisation of the process and plot it
out = gillespie(LV,10000)
op = par(mfrow=c(2,2))
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")

# use the "discretise" function to map it to an R "ts" object
plot(discretise(out,dt=0.01),plot.type="single",lty=c(1,2))
par(op)
}

\keyword{smfsb}
