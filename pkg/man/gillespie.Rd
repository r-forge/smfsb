\name{gillespie}
\alias{gillespie}
\title{Simulate a sample path from a stochastic kinetic model decrbribed by a stochastic Petri net}
\description{
This function simulates a single realisation from a discrete stochastic kinetic model described by a stochastic Petri net.
}
\usage{
gillespie(N, n, ...)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic Petri net. Should contain N$M, a vector representing the initial marking of the net, N$Pre, a matrix representing the LHS stoichiometries, N$Post, a matrix representing the RHS stoichiometries, and N$h, a function representing the rates of the reaction processes. N$h should have second argument t, a scalar representing the current simulation time (in the typical time-homogeneous case, N$h will ignore this argument), first argument x, a vector representing the current state of the system. N$h may posses additional arguments, representing reaction rates, for example.}
  \item{n}{An integer representing the number of events to simulate, excluding the initial state, N$M.}
  \item{...}{Additional arguments will be passed into the function N$h.}
}
\value{
An list with first component t, a vector of length n containing event times and second component x, a matrix with n+1 rows containing the state of the system. The ith row of x contains the state of the system prior to the ith event.
}
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}}, \code{\link{discretise}}, \code{\link{gillespied}}}
\examples{
# example

# Build a SPN for the Lotka-Volterra system
N=list()
N$M=c(x1=50,x2=100)
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(x,t,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# simulate a realisation of the process and plot it
out = gillespie(N,10000)
op = par(mfrow=c(2,2))
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")

# use the "discretise" function to map it to an R "ts" object
plot(discretise(out,dt=0.01),plot.type="single",lty=c(1,2))
par(op)
}

\keyword{smfsb}
