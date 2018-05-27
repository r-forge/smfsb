\name{StepCLE1D}
\alias{StepCLE1D}
\title{Create a function for advancing the state of an SPN by using a simple Euler-Maruyama discretisation of the CLE on a 1D regular grid}
\description{
This function creates a function for advancing the state of an SPN model
using a simple Euler-Maruyama discretisation of the CLE on a 1D regular grid. The resulting function (closure) can be
used in conjunction with other functions (such as \code{\link{simTs1D}})
for simulating realisations of SPN models in space and time.
}
\usage{
StepCLE1D(N,d,dt=0.01)
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
for a given molecule leaving the compartment is therefore twice this
value (as it can leave to the left or the right).}
\item{dt}{Time step for the Euler-Maruyama discretisation.}
}
\value{
An R function which can be used to advance the state of the SPN model
\code{N} by using a simple Euler-Maruyama algorithm. The function closure has
interface \code{function(x0,t0,deltat,...)}, where \code{x0} is a matrix
with rows corresponding to species and columns corresponding to voxels,
representing the initial condition, \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a matrix representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepGillespie1D}},\code{\link{StepCLE}},
  \code{\link{simTs1D}}, \code{\link{StepCLE2D}}}
\examples{
\donttest{
N=200
T=40
data(spnModels)
x0=matrix(0,nrow=2,ncol=N)
rownames(x0)=c("x1","x2")
x0[,round(N/2)]=LV$M
stepLV1D = StepCLE1D(LV,c(0.6,0.6),dt=0.05)
xx = simTs1D(x0,0,T,0.2,stepLV1D)
op=par(mfrow=c(1,2))
image(xx[1,,],main="Prey",xlab="Space",ylab="Time")
image(xx[2,,],main="Predator",xlab="Space",ylab="Time")
par(op)
}
}

\keyword{smfsb}
