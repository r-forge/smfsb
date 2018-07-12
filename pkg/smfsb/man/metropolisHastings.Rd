\name{metropolisHastings}
\alias{metropolisHastings}
\title{Run a Metropolis-Hastings MCMC algorithm for the parameters of a
  Bayesian posterior distribution}
\description{
Run a Metropolis-Hastings MCMC algorithm for the parameters of a
  Bayesian posterior distribution. Note that the algorithm carries over
  the old likelihood from the previous iteration, making it suitable for
  problems with expensive likelihoods, and also for "exact approximate"
  pseudo-marginal or particle marginal MH algorithms.
}
\usage{
metropolisHastings(init, logLik, rprop, dprop=function(new, old, ...){1},
                   dprior=function(x, ...){1}, iters=10000, thin=10,
                   verb=TRUE, debug=FALSE)
}
\arguments{
  \item{init}{An parameter vector with which to initialise the MCMC algorithm.}
  \item{logLik}{A function which takes a parameter (such as \code{init}) as
    its only required argument and returns the log-likelihood of the
    data. Note that it is fine for this to return the log of an unbiased
    estimate of the likelihood, in which case the algorithm will be an
    "exact approximate" pseudo-marginal MH algorithm.}
  \item{rprop}{A function which takes a parameter as its only required
    argument and returns a single sample from a proposal distribution.}
  \item{dprop}{A function which takes a new and old parameter as its
    first two required arguments and returns the (log) density of the
    new value conditional on the old. It should accept an optional
    parameter \code{log}, and must return the log-density when \code{log} is
    \code{TRUE}. Defaults to a flat function which causes this term to drop
    out of the acceptance probability. It is fine to use the default for
    _any_ _symmetric_ proposal, since the term will also drop out for
    any symmetric proposal.}
  \item{dprior}{A function which take a parameter as its only required
    argument and returns the (log) density of the parameter value under
    the prior. It should accept an optional
    parameter \code{log}, and must return the log-density when \code{log} is
    \code{TRUE}. Defaults to a flat function which causes this term to drop
    out of the acceptance probability. People often use a flat prior
    when they are trying to be "uninformative" or "objective", but this
    is slightly naive. In particular, what is "flat" is clearly
    dependent on the parametrisation of the model.}
  \item{iters}{The number of MCMC iterations required (_after_
    thinning).}
  \item{thin}{The required thinning factor. eg. only store every \code{thin}
    iterations.}
  \item{verb}{Boolean indicating whether some progress information
    should be printed to the console. Defaults to \code{TRUE}.}
  \item{debug}{Boolean indicating whether debugging information is required. Prints
    information about each iteration to console, to, eg., debug a crashing sampler.}
}
\value{
A matrix with rows representing samples from the posterior distribution.
}
\seealso{\code{\link{pfMLLik}}, \code{\link{StepGillespie}}, \code{\link{abcRun}},
 \code{\link{simTs}}, \code{\link{stepLVc}}, \code{\link{metrop}}}
\examples{
## First simulate some synthetic data
data = rnorm(250,5,2)
## Now use MH to recover the parameters
llik = function(x) { sum(dnorm(data,x[1],x[2],log=TRUE)) }
prop = function(x) { rnorm(2,x,0.1) }
prior = function(x, log=TRUE) {
    l = dnorm(x[1],0,100,log=TRUE) + dgamma(x[2],1,0.0001,log=TRUE)
    if (log) l else exp(l)
}
out = metropolisHastings(c(mu=1,sig=1), llik, prop,
                         dprior=prior, verb=FALSE)
out = out[1000:10000,]
mcmcSummary(out, truth=c(5,2), rows=2, plot=FALSE)
}

\keyword{smfsb}

