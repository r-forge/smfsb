\name{simpleEuler}
\alias{simpleEuler}
\title{Simulate a sample path from an ODE model}
\description{
This function integrates an Ordinary Differential Equation (ODE) model
using a simple first order Euler method. The function is pedagogic and
not intended for serious use. See the \code{\link[deSolve]{deSolve}}
package for better, more robust ODE solvers.
}
\usage{
simpleEuler(t=50, dt=0.001, fun, ic, ...)
}
\arguments{
  \item{t}{The length of the time interval over which the ODE model is to be integrated. Defaults to 50 time units.}
  \item{dt}{The step size to be used both for the time step of the Euler
    integration method and the recording interval for the output. It
    would probably be better to have separate parameters for these two
    things (see \code{\link{StepEuler}} and
    \code{\link{simTs}}). Defaults to 0.01 time units.} 
  \item{fun}{A vector-valued function representing the right hand side
    of the ODE model.
    The first argument is a vector representing the current state of the
    model, \code{x}. 
    The second argument of \code{fun} is the current simulation time,
    \code{t}. In the case of a homogeneous ODE model, this argument will
    be unused within the function. The function may have additional arguments, representing
    model parameters.
    The output of \code{fun} should be a vector of the same dimension as \code{x}.}
  \item{ic}{The initial conditions for the ODE model. This should be a vector of the same dimensions as the output from \code{fun}, and the second argument of \code{fun}.}
  \item{...}{Additional arguments will be passed into \code{fun}.}
}
\value{
An R \code{\link{ts}} object containing the sampled path of the model.
}
\seealso{\code{\link{rdiff}}, \code{\link{ts}}, \code{\link{StepEuler}},
\code{\link{simTs}}}
\examples{
# simple Lotka-Volterra example
lv <- function(x,t,k=c(k1=1,k2=0.1,k3=0.1))
{
        with(as.list(c(x,k)),{
                c( k1*x1 - k2*x1*x2 ,
                      k2*x1*x2 - k3*x2 )
        })
}
plot(simpleEuler(t=100,fun=lv,ic=c(x1=4,x2=10)),plot.type="single",lty=1:2)

# now an example which instead uses deSolve...
require(deSolve)
times = seq(0,50,by=0.01)
k = c(k1=1,k2=0.1,k3=0.1)
lvlist = function(t,x,k)
        list(lv(x,t,k))
plot(ode(y=c(x1=4,x2=10),times=times,func=lvlist,parms=k))

}

\keyword{smfsb}
