\name{StepODE}
\alias{StepODE}
\title{Create a function for advancing the state of an ODE model by using the deSolve package}
\description{
This function creates a function for advancing the state of an ODE model using an integration method from the \code{deSolve} package. The resulting function (closure) can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of ODE models.
}
\usage{
StepODE(RHSfun)
}
\arguments{
  \item{RHSfun}{A function representing the RHS of the ODE model. \code{RHSfun} should have prototype \code{RHSfun(t,x,parms,...)}, with \code{t} representing current system time, \code{x} representing current system state and \code{parm} representing the model parameters. The value of the function should be a vector of the same dimension as \code{x}, representing the infintesimal change in state.}
}
\value{
An R function which can be used to advance the state of the ODE model \code{RHSfun} by using an Euler method with step size \code{dt}. The function closure has interface \code{function(t0,x0,deltat,parms,...)}, where \code{t0} and \code{x0} represent the initial time and state, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEulerSPN}}, \code{\link{StepEuler}}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# Build a RHS for the Lotka-Volterra system
LVrhs <- function(t,x,parms)
{
        with(as.list(c(x,parms)),{
                c( c1*x1 - c2*x1*x2 ,
                      c2*x1*x2 - c3*x2 )
        })
}

# create a stepping function
stepLV=StepODE(LVrhs)

# step the function
print(stepLV(0,c(x1=50,x2=100),1,parms=c(c1=1,c2=0.005,c3=0.6)))

# integrate the process and plot it
out=simTs(0,c(x1=50,x2=100),50,0.1,stepLV,parms=c(c1=1,c2=0.005,c3=0.6))
plot(out,plot.type="single",lty=1:2)
}

\keyword{smfsb}
