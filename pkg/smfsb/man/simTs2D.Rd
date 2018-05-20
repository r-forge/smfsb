\name{simTs2D}
\alias{simTs2D}
\title{Simulate a model on a regular grid of times, using a function (closure) for advancing the state of the model}
\description{
This function simulates single realisation of a model on a 2D regular
spatial grid and regular grid of times using a function (closure) for advancing the state of the model, such as created by \code{\link{StepGillespie2D}}.
}
\usage{
simTs2D(x0,t0=0,tt=100,dt=0.1,stepFun,verb=FALSE,...)
}
\arguments{
  \item{x0}{The initial state of the process at time \code{t0}, a 3d array
    with dimensions corresponding to reacting species and two spatial dimensions.}
  \item{t0}{The initial time to be associated with the initial state \code{x0}. Defaults to 0.}
  \item{tt}{The terminal time of the simulation.}
  \item{dt}{The time step of the output. Note that this time step
    relates only to the recorded output, and has no bearing on the accuracy of the simulation process.}
  \item{stepFun}{A function (closure) for advancing the state of the
    process, such as produced by \code{\link{StepGillespie2D}}.}
  \item{verb}{Output progress to the console and graphics window (this function can be very slow).}
  \item{...}{Additional arguments will be passed to \code{stepFun}.}
}
\value{
An R 4d array representing the simulated process. The dimensions are
species, 2 space, and time.}
\seealso{\code{\link{StepGillespie2D}}, \code{\link{simTs1D}}}
\examples{
\donttest{
data(spnModels)
m=20; n=30; T=15
x0=array(0,c(2,m,n))
dimnames(x0)[[1]]=c("x1","x2")
x0[,round(m/2),round(n/2)]=LV$M
stepLV2D = StepGillespie2D(LV,c(0.6,0.6))
xx = simTs2D(x0,0,T,0.2,stepLV2D,verb=TRUE)
N = dim(xx)[4]
op=par(mfrow=c(1,2))
image(xx[1,,,N],main="Prey",xlab="Space",ylab="Time")
image(xx[2,,,N],main="Predator",xlab="Space",ylab="Time")
par(op)
}}

\keyword{smfsb}
