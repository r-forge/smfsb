\name{rfmc}
\alias{rfmc}
\title{Simulate a finite state space Markov chain}
\description{
This function simulates a single realisation from a discrete time Markov chain having a finite state space based on a given transition matrix.
}
\usage{
rfmc(n,P,pi0)
}
\arguments{
  \item{n}{The number of states to be sampled from the Markov chain, including the initial state, which will be sampled using \code{pi0}.}
  \item{P}{The transition matrix of the Markov chain. This is assumed to be a stochastic matrix, having non-negative elements and rows summing to one, though in fact, the rows will in any case be normalised by the sampling procedure.}
  \item{pi0}{A vector representing the probability distribution of the initial state of the Markov chain. If this vector is of length \code{r}, then the transition matrix \code{P} is assumed to be \code{r x r}. The elements of this vector are assumed to be non-negative and sum to one, though in fact, they will be normalised by the sampling procedure.}
}
\value{
An R \code{\link{ts}} object containing the sampled values from the Markov chain.
}
\seealso{\code{\link{rcfmc}}, \code{\link{ts}}}
\examples{
# example for sampling a finite Markov chain
P = matrix(c(0.9,0.1,0.2,0.8),ncol=2,byrow=TRUE)
pi0 = c(0.5,0.5)
samplepath = rfmc(200,P,pi0)
plot(samplepath)
summary(samplepath)
table(samplepath)
table(samplepath)/length(samplepath) # empirical distribution
# now compute the exact stationary distribution...
e = eigen(t(P))$vectors[,1]
e/sum(e)
}

\keyword{smfsb}
