\name{StepODE}
\alias{StepODE}
\title{Create a function for advancing the state of an ODE model by using the deSolve package}
\description{
This function creates a function for advancing the state of an ODE model
using an integration method from the \code{\link[deSolve]{deSolve}} package. The
resulting function (closure) can be used in conjunction with other
functions (such as \code{\link{simTs}}) for simulating realisations of
ODE models. This function is used similarly to \code{\link{StepEuler}},
but \code{\link{StepODE}} should be more accurate and efficient.
}
\usage{
StepODE(RHSfun)
}
\arguments{
  \item{RHSfun}{A function representing the RHS of the ODE model. \code{RHSfun} should have prototype \code{RHSfun(x,t,parms,...)}, with \code{t} representing current system time, \code{x} representing current system state and \code{parms} representing the model parameters. The value of the function should be a vector of the same dimension as \code{x}, representing the infinitesimal change in state.}
}
\value{
An R function which can be used to advance the state of the ODE model
\code{RHSfun} by using an efficient ODE solver. The function closure has interface \code{function(x0,t0,deltat,parms,...)}, where \code{t0} and \code{x0} represent the initial time and state, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepEuler}},
  \code{\link{simTs}}, \code{\link[deSolve]{ode}}}
\examples{
# Build a RHS for the Lotka-Volterra system
LVrhs <- function(x,t,parms)
{
        with(as.list(c(x,parms)),{
                c( c1*x1 - c2*x1*x2 ,
                      c2*x1*x2 - c3*x2 )
        })
}
# create a stepping function
stepLV = StepODE(LVrhs)
# step the function
print(stepLV(c(x1=50,x2=100),0,1,parms=c(c1=1,c2=0.005,c3=0.6)))
# integrate the process and plot it
out = simTs(c(x1=50,x2=100),0,50,0.1,stepLV,parms=c(c1=1,c2=0.005,c3=0.6))
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
