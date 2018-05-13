\name{rcfmc}
\alias{rcfmc}
\title{Simulate a continuous time finite state space Markov chain}
\description{
This function simulates a single realisation from a continuous time Markov chain having a finite state space based on a given transition rate matrix.
}
\usage{
rcfmc(n,Q,pi0)
}
\arguments{
  \item{n}{The number of states to be sampled from the Markov chain, including the initial state, which will be sampled using \code{pi0}.}
  \item{Q}{The transition rate matrix of the Markov chain, where each off-diagonal element \code{Q[i,j]} represents the rate of transition from state \code{i} to state \code{j}. This matrix is assumed to be square, having rows summing to zero.}
  \item{pi0}{A vector representing the probability distribution of the
    initial state of the Markov chain. If this vector is of length \code{r},
    then the transition matrix \code{P} is assumed to be \code{r x r}. The elements of this vector are assumed to be non-negative and sum to one, though in fact, they will be normalised by the sampling procedure.}
}
\value{
An R \code{\link{stepfun}} object containing the sampled path of the process.
}
\seealso{\code{\link{rfmc}}, \code{\link{stepfun}}}
\examples{
plot(rcfmc(20,matrix(c(-0.5,0.5,1,-1),ncol=2,byrow=TRUE),c(1,0)))
}

\keyword{smfsb}
