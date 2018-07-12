\name{abcRun}
\alias{abcRun}
\title{Run a set of simulations initialised with parameters sampled from a given prior distribution, and compute statistics required for an ABC analaysis}
\description{
Run a set of simulations initialised with parameters sampled from a given prior distribution, and compute statistics required for an ABC analaysis. Typically used to calculate "distances" of simulated synthetic data from observed data.
}
\usage{
abcRun(n, rprior, rdist)
}
\arguments{
  \item{n}{An integer representing the number of simulations to run.}
  \item{rprior}{A function without arguments generating a single parameter (vector) from prior distribution.}
  \item{rdist}{A function taking a parameter (vector) as argument and returning the required statistic of interest. This will typically be computed by first using the parameter to run a forward model, then computing required summary statistics, then computing a distance. See the example for details.}
}
\value{
A list with elements 'param' and 'dist'. These will be returned as matrices or vectors depending on whether the parameters and distances are scalars or vectors.
}
\seealso{\code{\link{pfMLLik}}, \code{\link{StepGillespie}}, \code{\link{abcSmc}},
 \code{\link{simTs}}, \code{\link{stepLVc}}}
\examples{
\donttest{
data(LVdata)
rprior <- function() { exp(c(runif(1, -3, 3),runif(1,-8,-2),runif(1,-4,2))) }
rmodel <- function(th) { simTs(c(50,100), 0, 30, 2, stepLVc, th) }
sumStats <- identity
ssd = sumStats(LVperfect)
distance <- function(s) {
    diff = s - ssd
    sqrt(sum(diff*diff))
}
rdist <- function(th) { distance(sumStats(rmodel(th))) }
out = abcRun(10000, rprior, rdist)
q=quantile(out$dist, c(0.01, 0.05, 0.1))
print(q)
accepted = out$param[out$dist < q[1],]
print(summary(accepted))
print(summary(log(accepted)))
}
}

\keyword{smfsb}
