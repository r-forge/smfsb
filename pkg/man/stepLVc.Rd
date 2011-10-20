\name{stepLVc}
\alias{stepLVc}
\title{A function for advancing the state of a Lotka-Volterra model by using the Gillespie algorithm}
\description{
A function for advancing the state of a Lotka-Volterra model by calling some C code implementing the Gillespie algorithm. The function can be used in conjunction with other functions (such as \code{\link{simTs}}) for simulating realisations of Lotka-Volterra models. Should be functionally identical to the function obtained by \code{data(spnModels)}, \code{stepLV=StepGillespie(LV)}, but much faster. 
}
\usage{
stepLVc(x0,t0,deltat,th=c(1,0.005,0.6))
}
\arguments{
\item{x0}{A vector representing the state of the system at the initial time, \code{t0}.}
\item{t0}{The time corresponding to the initial state, \code{x0}.}
\item{deltat}{The time in advance of the initial time at which the new state is required.}
\item{th}{A vector of length 3 representing the rate constants associated with the 3 LV reactions. Defaults to \code{c(1,0.005,0.6)}.}
}
\value{
A 2-vector representing the new state of the LV system.
}
\seealso{\code{\link{StepGillespie}}, \code{spnModels}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# load the LV model
data(spnModels)
# create a stepping function
stepLV = StepGillespie(LV)
# step the function
print(stepLV(c(x1=50,x2=100),0,1))
# simulate a realisation of the process and plot it
out = simTs(c(x1=50,x2=100),0,100,0.1,stepLV)
plot(out)
# now use "stepLVc" instead...
out = simTs(c(x1=50,x2=100),0,100,0.1,stepLVc)
plot(out)
}

\keyword{smfsb}
