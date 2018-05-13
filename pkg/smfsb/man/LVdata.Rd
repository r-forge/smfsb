\name{LVdata}
\alias{LVdata}
\alias{LVperfect}
\alias{LVprey}
\alias{LVnoise10Scale10}
\alias{LVpreyNoise10}
\alias{LVpreyNoise10Scale10}
\alias{LVnoise10}
\alias{LVnoise30}
\alias{LVnoise3010}
\alias{LVirregular}
\alias{LVirregularNoise10}
\docType{data}
\title{Example simulated time courses from a stochastic Lotka--Volterra model}
\description{
Collection of simulated time courses from a stochastic Lotka--Volterra
model.
\code{LVperfect} is direct output from a Gillespie simulation.
\code{LVprey} is the prey component.
\code{LVnoise10} has Gaussian noise with standard deviation 10 added.
\code{LVnoise30} has Gaussian noise with standard deviation 30 added.
\code{LVpreyNoise10} is the prey component with 10 SD noise added.
\code{LVnoise3010} has Gaussian noise added. The noise added to the prey
component has standard deviation 30 and the noise added to the predator
component has standard deviation 10.
\code{LVnoise10scale10} has Gaussian noise with standard deviation 10
added, and is then rescaled by a factor of 10 to mimic a scenario of an
uncalibrated measurement scale.
\code{LVirregular} is direct output from a Gillespie simulator, but on
an irregular time grid.
\code{LVirregularNoise10} is output on an irregular time grid with
Gaussian noise of standard deviation 10 added.
}
\usage{data(LVdata)}
\format{All datasets beginning
  \code{LVirregular} are R matrices such as output by
  \code{\link{simTimes}}, and the rest are R \code{\link{ts}} objects
  such as output by \code{\link{simTs}}.} 
\keyword{smfsb}
\keyword{data}
\keyword{datasets}

