\name{pfMLLik1}
\alias{pfMLLik1}
\title{Create a function for computing the log of an unbiased estimate of marginal likelihood of a time course data set}
\description{
Create a function for computing the log of an unbiased estimate of marginal likelihood of a time course data set using a simple bootstrap particle filter. This version does not use the "log-sum-exp trick" for avoiding numerical underflow. 
See \code{\link{pfMLLik}} for a version which does.
}
\usage{
pfMLLik1(n,simx0,t0,stepFun,dataLik,data)
}
\arguments{
  \item{n}{An integer representing the number of particles to use in the particle filter.}
  \item{simx0}{A function with interface \code{simx0(n,t0,...)}, where \code{n} is the number of rows of a matrix and \code{t0} is a time at which to simulate from an initial distribution for the state of the particle filter. The return value should be a matrix whose rows are random samples from this distribution. The function therefore represents a prior distribution on the initial state of the Markov process.}
  \item{t0}{The time corresponding to the starting point of the Markov process. Can be no bigger than the smallest observation time.}
  \item{stepFun}{A function for advancing the state of the Markov process, such as returned by \code{\link{StepGillespie}}.}
  \item{dataLik}{A function with interface
    \code{dataLik(x,t,y,log=TRUE,...)}, where \code{x} and \code{t}
    represent the true state and time of the process, and \code{y} is
    the observed data. The return value should be the (log of the)
    likelihood of the observation. The function therefore represents the
  observation model.}
  \item{data}{A timed data matrix representing the observations, such as produced by \code{\link{simTimes}} or \code{\link{as.timedData}}.}
}
\value{
An R function with interface \code{(...)} which evaluates to the log of the particle filter's unbiased estimate of the marginal likelihood of the data.
}
\seealso{\code{\link{pfMLLik}}, \code{\link{StepGillespie}}, \code{\link{as.timedData}},
  \code{\link{simTimes}}, \code{\link{stepLVc}}}
\examples{
noiseSD=5
# first simulate some data
truth=simTs(c(x1=50,x2=100),0,20,2,stepLVc)
data=truth+rnorm(prod(dim(truth)),0,noiseSD)
data=as.timedData(data)
# measurement error model
dataLik <- function(x,t,y,log=TRUE,...)
{
    ll=sum(dnorm(y,x,noiseSD,log=TRUE))
    if (log)
        return(ll)
    else
        return(exp(ll))
}
# now define a sampler for the prior on the initial state
simx0 <- function(N,t0,...)
{
    mat=cbind(rpois(N,50),rpois(N,100))
    colnames(mat)=c("x1","x2")
    mat
}
mLLik=pfMLLik1(1000,simx0,0,stepLVc,dataLik,data)
print(mLLik())
print(mLLik(th=c(th1 = 1, th2 = 0.005, th3 = 0.6)))
print(mLLik(th=c(th1 = 1, th2 = 0.005, th3 = 0.5)))
}

\keyword{smfsb}
