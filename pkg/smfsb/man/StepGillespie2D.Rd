\name{StepGillespie2D}
\alias{StepGillespie2D}
\title{Create a function for advancing the state of an SPN by using the
  Gillespie algorithm on a 2D regular grid}
\description{
This function creates a function for advancing the state of an SPN model
using the Gillespie algorithm. The resulting function (closure) can be
used in conjunction with other functions (such as \code{\link{simTs2D}})
for simulating realisations of SPN models in space and time.
}
\usage{
StepGillespie2D(N,d)
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
    \code{N$h} may possess additional arguments, representing reaction
    rates, for example. \code{N} does not need to contain an initial
    marking, \code{N$M}. \code{N$M} will be ignored by most functions
    which use the resulting function closure.}
  \item{d}{A vector of diffusion coefficients - one coefficient for each
  reacting species, in order. The coefficient is the reaction rate for a
reaction for a molecule moving into an adjacent compartment. The hazard
for a given molecule leaving the compartment is therefore four times this
value (as it can leave in one of 4 directions).}
}
\value{
An R function which can be used to advance the state of the SPN model
\code{N} by using the Gillespie algorithm. The function closure has
interface \code{function(x0,t0,deltat,...)}, where \code{x0} is a 3d array
with dimensions corresponding to species followed by two spatial dimensions,
representing the initial condition, \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns an array representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepGillespie}},
  \code{\link{simTs2D}}, \code{\link{StepGillespie1D}}}
\examples{
\donttest{
data(spnModels)
m=20; n=30; T=10
x0=array(0,c(2,m,n))
dimnames(x0)[[1]]=c("x1","x2")
x0[,round(m/2),round(n/2)]=LV$M
stepLV2D = StepGillespie2D(LV,c(0.6,0.6))
xx = simTs2D(x0,0,T,0.2,stepLV2D,verb=TRUE)
N = dim(xx)[4]
op=par(mfrow=c(1,2))
image(xx[1,,,N],main="Prey",xlab="Space",ylab="Time")
image(xx[2,,,N],main="Predator",xlab="Space",ylab="Time")
par(op)
}}

\keyword{smfsb}
