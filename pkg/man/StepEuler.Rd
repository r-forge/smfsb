\name{StepEuler}
\alias{StepEuler}
\title{Create a function for advancing the state of an ODE model by using a simple Euler integration method}
\description{
This function creates a function for advancing the state of an ODE model
using a simple Euler integration method. The resulting function
(closure) can be used in conjunction with other functions (such as
\code{\link{simTs}}) for simulating realisations of ODE models. This
function is intended to be pedagogic. See \code{\link{StepODE}} for a
more accurate integration function.
}
\usage{
StepEuler(RHSfun,dt=0.01)
}
\arguments{
  \item{RHSfun}{A function representing the RHS of the ODE
    model. \code{RHSfun} should have prototype \code{RHSfun(x,t,...)},
    with  \code{x} representing current system state and \code{t} representing current system time. The value of the function should be a vector of the same dimension as \code{x}, representing the infinitesimal change in state.}
  \item{dt}{Time step to be used by the simple Euler integration method. Defaults to 0.01.}
}
\value{
An R function which can be used to advance the state of the ODE model \code{RHSfun} by using an Euler method with step size \code{dt}. The function closure has interface \code{function(x0,t0,deltat,...)}, where \code{t0} and \code{x0} represent the initial time and state, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepODE}}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# Build a RHS for the Lotka-Volterra system
LVrhs <- function(x,t,th=c(c1=1,c2=0.005,c3=0.6))
{
        with(as.list(c(x,th)),{
                c( c1*x1 - c2*x1*x2 ,
                      c2*x1*x2 - c3*x2 )
        })
}
# create a stepping function
stepLV = StepEuler(LVrhs)
# step the function
print(stepLV(c(x1=50,x2=100),0,1))
# integrate the process and plot it
out = simTs(c(x1=50,x2=100),0,20,0.1,stepLV)
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
