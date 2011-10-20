\name{normgibbs}
\alias{normgibbs}
\title{A simple Gibbs sampler for Bayesian inference for the mean and precision of a normal random sample}
\description{
This function runs a simple Gibbs sampler for the Bayesian posterior distribution of the mean and precision given a normal random sample.
}
\usage{
normgibbs(N, n, a, b, cc, d, xbar, ssquared)
}
\arguments{
  \item{N}{The number of iterations of the Gibbs sampler.}
  \item{n}{The sample size of the normal random sample.}
  \item{a}{The shape parameter of the gamma prior on the sample precision.}
  \item{b}{The scale parameter of the gamma prior on the sample precision.}
  \item{cc}{The mean of the normal prior on the sample mean.}
  \item{d}{The precision of the normal prior on the sample mean.}
  \item{xbar}{The sample mean of the data. eg. \code{mean(x)} for a vector \code{x}.}
  \item{ssquared}{The sum of squares of the data. eg. \code{sum(x*x)} for a vector \code{x}.}
}
\value{
An R matrix object containing the samples of the Gibbs sampler.
}
\seealso{\code{\link{rcfmc}}, \code{\link{metrop}}, \code{\link{mcmcSummary}}}
\examples{
postmat=normgibbs(N=1100,n=15,a=3,b=11,cc=10,d=1/100,xbar=25,ssquared=20)
postmat=postmat[101:1100,]
op=par(mfrow=c(3,3))
plot(postmat)
plot(postmat,type="l")
plot.new()
plot(ts(postmat[,1]))
plot(ts(postmat[,2]))
plot(ts(sqrt(1/postmat[,2])))
hist(postmat[,1],30)
hist(postmat[,2],30)
hist(sqrt(1/postmat[,2]),30)
par(op)
}

\keyword{smfsb}
