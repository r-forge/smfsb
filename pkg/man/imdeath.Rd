\name{imdeath}
\alias{imdeath}
\title{Simulate a sample path from the homogeneous immigration-death process}
\description{
This function simulates a single realisation from a time-homogeneous immigration-death process.
}
\usage{
imdeath(n=20,x0=0,lambda=1,mu=0.1)
}
\arguments{
  \item{n}{The number of states to be sampled from the process, not including the initial state, \code{x0}}
  \item{x0}{The initial state of the process, which defaults to zero.}
  \item{lambda}{The rate at which new individual immigrate into the population. Defaults to 1.}
  \item{mu}{The rate at which individuals within the population die, independently of all other individuals. Defaults to 0.1.}
}
\value{
An R \code{\link{stepfun}} object containing the sampled path of the process.
}
\seealso{\code{\link{rcfmc}}, \code{\link{rdiff}},
  \code{\link{stepfun}}, \code{\link{gillespie}}}
\examples{
plot(imdeath(50))
}

\keyword{smfsb}
