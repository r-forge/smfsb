\name{gillespied}
\alias{gillespied}
\title{Simulate a sample path from a stochastic kinetic model described by a stochastic Petri net}
\description{
This function simulates a single realisation from a discrete stochastic kinetic model described by a stochastic Petri net and discretises the output onto a regular time grid.
}
\usage{
gillespied(N, T=100, dt=1, ...)
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
    \code{N$h} may posses additional arguments, representing reaction rates, for example.
}
  \item{T}{The required length of simulation time. Defaults to 100 time units.}
  \item{dt}{The grid size for the output. Note that this parameter simply determines the volume of output. It has no bearing on the correctness of the simulation algorithm. Defaults to one time unit.}
  \item{...}{Additional arguments will be passed into the function \code{N$h}.}
}
\value{
An R \code{\link{ts}} object containing the simulated realisation of the process. }
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}},
  \code{\link{discretise}}, \code{\link{gillespie}}, \code{\link{StepGillespie}}}
\examples{
# load LV model
data(spnModels)
# simulate and plot a realisation
plot(gillespied(LV,T=100,dt=0.01))
}

\keyword{smfsb}
