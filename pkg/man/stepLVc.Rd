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
\item{x0}{A vector representing the state of the system at the inital time, \code{tt}.}
\item{t0}{The time corresponding to the initial state, \code{x}.}
\item{deltat}{The time in advance of the initial time at which the new state is required.}
\item{th}{A vector of length 3 representing the rate constants associated with the 3 LV reactions. Defaults to \code{c(1,0.005,0.6)}.}
}
\value{
A 2-vector representing the new state of the LV system.
}
\seealso{\code{\link{StepGillespie}}, \code{data(spnModels)}, \code{\link{simTs}}, \code{\link{simSample}}}
\examples{
# example
# Build a SPN for the Lotka-Volterra system
N=list()
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(x,t,th=c(th1=1,th2=0.005,th3=0.6))
{
 with(as.list(c(x,th)),{
         return(c(th1*x1, th2*x1*x2, th3*x2 ))
        })
}

# create a stepping function
stepLV = StepGillespie(N)

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
