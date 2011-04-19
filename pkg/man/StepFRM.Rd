\name{StepFRM}
\alias{StepFRM}
\title{Create a function for advancing the state of an SPN by using Gillespie's first reaction method}
\description{
This function creates a function for advancing the state of an SPN model using Gillespie's first reaction method. The resulting function (closure) can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of SPN models.
}
\usage{
StepFRM(N)
}
\arguments{
  \item{N}{An R list with named components representing a stochastic Petri net. Should contain \code{N$Pre}, a matrix representing the LHS stoichiometries, \code{N$Post}, a matrix representing the RHS stoichiometries, and \code{N$h}, a function representing the rates of the reaction processes. \code{N$h} should have first argument \code{t}, a scalar representing the current simulation time (in the typical time-homogeneous case, \code{N$h} will ignore this argument), second argument \code{x}, a vector representing the current state of the system. \code{N$h} may posses additional arguments, representing reaction rates, for example. \code{N} does not need to contain an initial marking, \code{N$M}. \code{N$M} will be ignored by most functions which use the resulting function closure.}
}
\value{
An R function which can be used to advance the state of the SPN model \code{N} by using Gillespie's first reaction method. The function closure has interface \code{function(x0,t0,delat,...)}, where \code{x0} and \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{gillespied}}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# example
# Build a SPN for the Lotka-Volterra system
N=list()
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(t,x,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# create a stepping function
stepLV=StepFRM(N)

# step the function
print(stepLV(0,c(x1=50,x2=100),1))

# simulate a realisation of the process and plot it
out=simTs(0,c(x1=50,x2=100),100,0.1,stepLV)
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
