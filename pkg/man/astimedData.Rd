\name{as.timedData}
\alias{as.timedData}
\title{Convert a time series object to a timed data matrix}
\description{
This function converts a time series object to a timed data matrix,
similar to that produced by \code{\link{simTimes}}. The main purpose is
for passing data to the function \code{\link{pfMLLik}}, which expects
data encoded in this format.
}
\usage{
as.timedData(timeseries)
}
\arguments{
  \item{timeseries}{An R timeseries object, such as produced by the functions \code{\link{ts}} or \code{\link{simTs}}.}
}
\value{
An R matrix object with row names corresponding to observation times, similar to that produced by \code{\link{simTimes}}.
}
\seealso{\code{\link{simTs}}, \code{\link{ts}}, \code{\link{simTimes}}, \code{\link{pfMLLik}}}
\examples{
truth=simTs(c(x1=50,x2=100),0,20,2,stepLVc)
simData=truth+rnorm(prod(dim(truth)),0,5)
timedData=as.timedData(simData)
print(timedData)
}

\keyword{smfsb}
