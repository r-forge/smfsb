\name{metrop}
\alias{metrop}
\title{Run a simple Metropolis sampler with standard normal target and uniform innovations}
\description{
This function runs a simple Metropolis sampler with standard normal target distribution and uniform innovations.
}
\usage{
metrop(n, alpha) 
}
\arguments{
  \item{n}{The number of iterations of the Metropolis sampler.}
  \item{alpha}{The tuning parameter of the sampler. The innovations of the sampler are of the form U(-alpha,alpha).}
}
\value{
An R vector containing the output of the sampler.
}
\seealso{\code{\link{normgibbs}}}
\examples{
normvec=metrop(1000,1)
op=par(mfrow=c(2,1))
plot(ts(normvec))
hist(normvec,20)
par(op)
}

\keyword{smfsb}
