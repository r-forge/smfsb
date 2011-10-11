\name{discretise}
\alias{discretise}
\title{Discretise output from a discrete event simulation algorithm}
\description{
This function discretise output from a discrete event simulation algorithm such as \code{\link{gillespie}} onto a regular time grid, and returns the results as an R \code{\link{ts}} object.
}
\usage{
discretise(out, dt=1, start=0)
}
\arguments{
  \item{out}{A list containing discrete event simulation output in the form of that produced by \code{\link{gillespie}}.}
  \item{dt}{The time step required for the output of the discretisation process. Defaults to one time unit.}
  \item{start}{The start time for the output. Defaults to zero.}
}
\value{
An R \code{\link{ts}} object containing the discretised output.
}
\seealso{\code{\link{simpleEuler}}, \code{\link{rdiff}}, \code{\link{gillespie}}, \code{\link{gillespied}}, \code{\link{ts}}}
\examples{
# load LV model
data(spnModels)
# simulate a realisation of the process and plot it
out = gillespie(LV,10000)
op=par(mfrow=c(2,2))
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")

# use the "discretise" function to map it to an R "ts" object
plot(discretise(out,dt=0.01),plot.type="single",lty=c(1,2))
par(op)
}

\keyword{smfsb}
