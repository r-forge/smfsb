\name{simTs1D}
\alias{simTs1D}
\title{Simulate a model on a regular grid of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates single realisation of a model on a 1D regular
spatial grid and regular grid of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie1D}}.
}
\usage{
simTs1D(x0,t0=0,tt=100,dt=0.1,stepFun,verb=FALSE,...)
}
\arguments{
  \item{x0}{The initial state of the process at time \code{t0}, a matrix
    with rows corresponding to reacting species and columns
    corresponding to spatial location.}
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{tt}{The terminal time of the simulation.}
  \item{dt}{The time step of the output. Note that this time step
    relates only to the recorded output, and has no bearing on the accuracy of the simulation process.}
  \item{stepFun}{A function (closure) for advancing the state of the
    process, such as produced by \code{\link{StepGillespie1D}}.}
  \item{verb}{Output progress to the console (this function can be very slow).}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R 3d array representing the simulated process. The dimensions are
species, space, and time.}
\seealso{\code{\link{StepGillespie1D}}, \code{\link{simTs}}}
\examples{
\donttest{
data(spnModels)
N=20; T=30
x0=matrix(0,nrow=2,ncol=N)
rownames(x0)=c("x1","x2")
x0[,round(N/2)]=LV$M
stepLV1D = StepGillespie1D(LV,c(0.6,0.6))
xx = simTs1D(x0,0,T,0.2,stepLV1D,verb=TRUE)
op=par(mfrow=c(1,2))
image(xx[1,,],main="Prey",xlab="Space",ylab="Time")
image(xx[2,,],main="Predator",xlab="Space",ylab="Time")
par(op)
}
}

\keyword{smfsb}
