\name{StepCLE2D}
\alias{StepCLE2D}
\title{Create a function for advancing the state of an SPN by using a simple Euler-Maruyama discretisation of the CLE on a 2D regular grid}
\description{
This function creates a function for advancing the state of an SPN model
using a simple Euler-Maruyama discretisation of the CLE on a 2D regular grid. The resulting function (closure) can be
used in conjunction with other functions (such as \code{\link{simTs2D}})
for simulating realisations of SPN models in space and time.
}
\usage{
StepCLE2D(N,d,dt=0.01)
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
\item{dt}{Time step for the Euler-Maruyama discretisation.}
}
\value{
An R function which can be used to advance the state of the SPN model
\code{N} by using a simple Euler-Maruyama algorithm. The function closure has
interface \code{function(x0,t0,deltat,...)}, where \code{x0} is a 3D array
with rows corresponding to species and columns corresponding to voxels,
representing the initial condition (with dimensions species, x, and y), \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a matrix representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepGillespie2D}},\code{\link{StepCLE}},
  \code{\link{simTs1D}}, \code{\link{StepCLE1D}}}
\examples{
\donttest{
m=150
n=100
T=15
data(spnModels)
x0=array(0,c(2,m,n))
dimnames(x0)[[1]]=c("x1","x2")
x0[,round(m/2),round(n/2)]=LV$M
stepLV2D = StepCLE2D(LV,c(0.6,0.6),dt=0.05)
xx = simTs2D(x0,0,T,0.5,stepLV2D,verb=TRUE)
N = dim(xx)[4]
op=par(mfrow=c(1,2))
image(xx[1,,,N],main="Prey",xlab="Space",ylab="Time")
image(xx[2,,,N],main="Predator",xlab="Space",ylab="Time")
par(op)
}
}

\keyword{smfsb}
