\name{gillespied}
\alias{gillespied}
\title{Simulate a sample path from a stochastic kinetic model decrbribed by a stochastic Petri net}
\description{
This function simulates a single realisation from a discrete stochastic kinetic model described by a stochastic Petri net and discretises the output onto a regular time grid.
}
\usage{
gillespied(N, T=100, dt=1, ...)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic Petri net. Should contain N$M, a vector representing the initial marking of the net, N$Pre, a matrix representing the LHS stoichiometries, N$Post, a matrix representing the RHS stoichiometries, and N$h, a function representing the rates of the reaction processes. N$h should have first argument t, a scalar representing the current simulation time (in the typical time-homogeneous case, N$h will ignore this argument), second argument x, a vector representing the current state of the system. N$h may posses additional arguments, representing reaction rates, for example.}
  \item{T}{The required length of simulation time. Defaults to 100 time units.}
  \item{dt}{The grid size for the output. Note that this parameter simply determines the volume of output. It has no bearing on the correctness of the simulation algorithm. Defaults to one time unit.}
  \item{...}{Additional arguments will be passed into the function N$h.}
}
\value{
An R \code{\link{ts}} object containing the simulated realisation of the process. }
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}}, \code{\link{discretise}}, \code{\link{gillespie}}}
\examples{
# example

N=list()
N$M=c(x1=50,x2=100)
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(t,x,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
  return(c(th1*x1, th2*x1*x2, th3*x2 ))
 })
}

plot(gillespied(N,T=100,dt=0.01))
}

\keyword{smfsb}
