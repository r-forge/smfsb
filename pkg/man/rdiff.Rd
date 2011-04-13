\name{rdiff}
\alias{rdiff}
\title{Simulate a sample path from a univariate diffusion process}
\description{
This function simulates a single realisation from a time-homogeneous univariate diffusion process.
}
\usage{
rdiff(afun, bfun, x0 = 0, t = 50, dt = 0.01, ...)
}
\arguments{
  \item{afun}{A scalar-valued function representing the infintesimal mean of the diffusion process. The first argument of afun is the current state of the process.}
  \item{bfun}{A scalar-valued function representing the infintesimal standard deviation of the process. The first argument of bfun is the current state of the process.}
  \item{x0}{The initial state of the diffusion process.}
  \item{t}{The length of the time interval over which the diffusion process is to be simulated. Defaults to 50 time units.}
  \item{dt}{The step size to be used both for the time step of the Euler integration method and the recording interval for the output. It would probably be better to have separate parameters for these two things. Defaults to 0.01 time units.}
  \item{...}{Additional arguments will be passed into afun and bfun.}
}
\value{
An R ts object containing the sampled path of the process.
}
\seealso{\code{\link{rcfmc}}}
\examples{
# example
# simulate a diffusion approximation to an immigration-death process
# infintesimal mean
afun<-function(x,lambda,mu)
{
        lambda-mu*x
}
# infintesimal standard deviation
bfun<-function(x,lambda,mu)
{
        sqrt(lambda+mu*x)
}
# plot a sample path
plot(rdiff(afun,bfun,lambda=1,mu=0.1,t=30))
}

\keyword{smfsb}
