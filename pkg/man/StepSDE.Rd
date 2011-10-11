\name{StepSDE}
\alias{StepSDE}
\title{Create a function for advancing the state of an SDE model by using a simple Euler-Maruyama integration method}
\description{
This function creates a function for advancing the state of an SDE model using a simple Euler-Maruyama integration method. The resulting function (closure) can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of SDE models.
}
\usage{
StepSDE(drift,diffusion,dt=0.01)
}
\arguments{
  \item{drift}{A function representing the drift vector of the SDE model
    (corresponding roughly to the RHS of an ODE model). \code{drift}
    should have prototype \code{drift(x,t,...)}, with
    \code{x} representing current system state and
    \code{t} representing current system time.
    The value of the function should be a vector of the same dimension as \code{x}, representing the infinitesimal mean of the Ito SDE.}
  \item{diffusion}{A function representing the diffusion matrix of the SDE model (the square root of the infinitesimal variance matrix). \code{diffusion} should have prototype \code{diffusion(x,t,...)}, with \code{x} representing current system state and \code{t} representing current system time. The value of the function should be a square matrix with both dimensions the same as the length of \code{x}.}
  \item{dt}{Time step to be used by the simple Euler-Maruyama integration method. Defaults to 0.01.}
}
\value{
An R function which can be used to advance the state of the SDE model with given drift vector and diffusion matrix, by using an Euler-Maruyama method with step size \code{dt}. The function closure has interface \code{function(x0,t0,deltat,...)}, where \code{x0} and \code{t0} represent the initial state and time, and \code{deltat} represents the amount of time by which the process should be advanced. The function closure returns a vector representing the simulated state of the system at the new time.}
\seealso{\code{\link{StepEuler}}, \code{\link{StepCLE}}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# Immigration-death diffusion approx with death rate a CIR process
myDrift <- function(x,t,th=c(lambda=1,alpha=1,mu=0.1,sigma=0.1))
     {
             with(as.list(c(x,th)),{
                     c( lambda - x*y ,
                           alpha*(mu-y) )
             })
     }
myDiffusion <- function(x,t,th=c(lambda=1,alpha=1,mu=0.1,sigma=0.1))
     {
             with(as.list(c(x,th)),{
                     matrix(c( sqrt(lambda + x*y) , 0,
                           0, sigma*sqrt(y) ),ncol=2,nrow=2,byrow=TRUE)
             })
     }
# create a stepping function
stepProc = StepSDE(myDrift,myDiffusion)
# integrate the process and plot it
out = simTs(c(x=5,y=0.1),0,20,0.1,stepProc)
plot(out)
}

\keyword{smfsb}
